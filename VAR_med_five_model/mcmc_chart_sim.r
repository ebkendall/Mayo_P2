library(matrixStats)
library(plotrix)

simulation = T

data_num = 4

trialNum = 1
sampNum = 4
itNum = 1
all_seeds = F
long_chain = F

if(all_seeds) {
    seed_list = 1:3
} else {
    seed_list = 1
}

# Load the model output -------------------------------------------------------
B_chain   = NULL
Hr_chain  = NULL
Map_chain = NULL
Hc_chain  = NULL
La_chain  = NULL

for(seed_num in 1:length(seed_list)) {
    print(seed_num)
    
    if(long_chain) {
        it_seq = 1:itNum
    } else {
        it_seq = itNum
    }
    
    for(it in it_seq) {
        file_name = paste0('Model_out/mcmc_out_interm_', seed_num,'_', trialNum,'it', it, '_samp', sampNum, '_sim.rda')
        load(file_name)
        print(file_name)
        
        if(it == 1) {
            B_chain   = mcmc_out_temp$B_chain[500:1000, ]
            Hr_chain  = mcmc_out_temp$hr_chain[500:1000, ]
            Map_chain = mcmc_out_temp$bp_chain[500:1000, ]
            Hc_chain  = mcmc_out_temp$hc_chain[500:1000, ]
            La_chain  = mcmc_out_temp$la_chain[500:1000, ]
        } else {
            B_chain   = rbind(B_chain, mcmc_out_temp$B_chain)
            Hr_chain  = mcmc_out_temp$hr_chain
            Map_chain = mcmc_out_temp$bp_chain
            Hc_chain  = mcmc_out_temp$hc_chain
            La_chain  = mcmc_out_temp$la_chain
        }
        rm(mcmc_out_temp)
    }
}
#  ----------------------------------------------------------------------------

# Mode of the state sequences -------------------------------------------------
Mode <- function(x) {
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
}

state_seq_mode = apply(B_chain, 2, Mode)

load(paste0('Data_sim/use_data1_', data_num, '.rda'))
load(paste0('Data_sim/Dn_omega_sim_', data_num, '.rda'))
load(paste0('Data_sim/omega_i_mat_', data_num, '.rda'))
Dn_omega = Dn_omega_sim
rm(Dn_omega_sim)
EIDs = unique(use_data[,'EID'])

ss_truth = use_data[,"b_true"]
load('Data_updates/hr_map_names1.rda')

# ------------------------------------------------------------------------------
# Function to change transparency of colors # ----------------------------------
# ------------------------------------------------------------------------------
makeTransparent = function(..., alpha=0.35) {
    
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    
    alpha = floor(255*alpha)  
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    
    .makeTransparent = function(col, alpha) {
        rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    
    return(newColor)
    
}

eid_poor = NULL

print("Summary of identifying correct states with mode")
if(length(ss_truth) != length(state_seq_mode)) {
    print("ERROR")
} else {
    print(sum(ss_truth == state_seq_mode) / length(state_seq_mode))   
}

print("Average (across subject) of proportion of correct states with mode")

prop_sub = rep(0, length(EIDs))
for(j in 1:length(EIDs)) {
    sub_ind_j = which(use_data[,"EID"] == EIDs[j])
    ss_truth_j = ss_truth[sub_ind_j]
    state_seq_mode_j = state_seq_mode[sub_ind_j]
    
    prop_sub[j] = sum(ss_truth_j == state_seq_mode_j) / length(state_seq_mode_j)
}

print(summary(prop_sub))
print(sort(prop_sub))

eid_poor = EIDs[prop_sub < 0.9]

state_2_info = which(ss_truth == 2)
state_2_truth = ss_truth[state_2_info]
state_2_mode  = state_seq_mode[state_2_info]

print("Correct identification of bleeding")
print(sum(state_2_truth == state_2_mode) / length(state_2_mode))  

# Sensitivity of state 2: Pr(predict S2 | true S2)
predict_at_true_S2 = state_seq_mode[ss_truth == 2]
print(paste0("Sensitivity of S2 = ", mean(predict_at_true_S2 == 2)))

# Specificity of state 2: Pr(predict not S2 | true not S2)
predict_not_S2 = state_seq_mode[ss_truth != 2]
print(paste0("Specificity of S2 = ", mean(predict_not_S2 != 2)))

# ------------------------------------------------------------------------------ 
# Model evaluation plots -------------------------------------------------------
# ------------------------------------------------------------------------------
pdf_title = NULL
if(all_seeds) {
    pdf_title = paste0('Plots/model_eval_', trialNum, '_it', itNum, '_samp', sampNum, '_sim.pdf') 
} else {
    if(long_chain) {
        pdf_title = paste0('Plots/model_eval_', trialNum, '_it', itNum, 
                           '_seed',seed_list, '_samp', sampNum, '_sim_long.pdf')
    } else {
        pdf_title = paste0('Plots/model_eval_', trialNum, '_it', itNum, 
                           '_seed',seed_list, '_samp', sampNum, '_sim.pdf')
    }
}
pdf(pdf_title)
panel_dim = c(4,1)
inset_dim = c(0,-.18)
par(mfrow=panel_dim, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
    print(which(EIDs == i))
    indices_i = (use_data[,'EID']==i)
    n_i = sum(indices_i)
    
    t_grid = round(use_data[indices_i, 'time'] / 60, digits = 3)
    t_grid_bar = 1:length(t_grid)
    rbc_times_bar = which(use_data[use_data[,'EID']==i, 'RBC_ordered'] != 0)
    if(simulation) {
        rbc_admin = c(head(use_data[use_data[,'EID']==i, "n_RBC_admin"], 1),
                      diff(use_data[use_data[,'EID']==i, "n_RBC_admin"]))
        rbc_admin_times_bar = which(rbc_admin != 0)
    } else {
        rbc_admin_times_bar = which(use_data[use_data[,'EID']==i, 'RBC_admin'] != 0)   
    }
    rbc_times = t_grid[rbc_times_bar]
    rbc_admin_times = t_grid[rbc_admin_times_bar]
    
    # Put this on the correct scale as the t_grid
    b_i = use_data[ indices_i,'b_true']
    to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
    to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
    to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
    to_s4 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==4]
    to_s5 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==5]
    
    if(b_i[1] == 1) {
        to_s1 = c(to_s1, 1)
    } else if(b_i[1] == 2) {
        to_s2 = c(to_s2, 1)
    } else if(b_i[1] == 3) {
        to_s3 = c(to_s3, 1)
    } else if(b_i[1] == 4) {
        to_s4 = c(to_s4, 1)
    } else {
        to_s5 = c(to_s5, 1)
    }
    
    if(length(unique(b_i)) > 1) {
        if(length(to_s1) > 0) {
            rect_coords = data.frame(s = 1, t = to_s1)
        }
        
        if(length(to_s2) > 0) {
            s2_coords = data.frame(s = 2, t = to_s2)
            if(length(to_s1) > 0) {
                rect_coords = rbind(rect_coords, s2_coords)
            } else {
                rect_coords = s2_coords
            }
        }
        
        if(length(to_s3) > 0) {
            s3_coords = data.frame(s = 3, t = to_s3)
            if(length(to_s1) > 0 || length(to_s2) > 0) {
                rect_coords = rbind(rect_coords, s3_coords)
            } else {
                rect_coords = s3_coords
            }
        }
        
        if(length(to_s4) > 0) {
            s4_coords = data.frame(s = 4, t = to_s4)
            if(length(to_s1) > 0 || length(to_s2) > 0 || length(to_s3) > 0) {
                rect_coords = rbind(rect_coords, s4_coords)
            } else {
                rect_coords = s4_coords
            }
        }
        
        if(length(to_s5) > 0) {
            s5_coords = data.frame(s = 5, t = to_s5)
            if(length(to_s1) > 0 || length(to_s2) > 0 || 
               length(to_s3) > 0 || length(to_s4) > 0) {
                rect_coords = rbind(rect_coords, s5_coords)
            } else {
                rect_coords = s5_coords
            }
        }
        
        if(!(n_i %in% rect_coords$t)) rect_coords = rbind(rect_coords, c(b_i[n_i], n_i))
        # Add one row for visuals
        rect_coords = rbind(rect_coords, c(b_i[n_i], n_i+1))
        rect_coords$t = rect_coords$t - 1
        rect_coords = rect_coords[order(rect_coords$t), ]
        col_vec = c('dodgerblue', 'firebrick1', 'yellow2', 
                    'green', 'darkgray')[rect_coords$s]
        col_vec = makeTransparent(col_vec, alpha = 0.35)   
    } else {
        rect_coords = data.frame(s = rep(b_i[1], 2), t = c(1,n_i+1))
        rect_coords$t = rect_coords$t - 1
        rect_coords = rect_coords[order(rect_coords$t), ]
        col_vec = c('dodgerblue', 'firebrick1', 'yellow2', 
                    'green', 'darkgray')[rect_coords$s]
        col_vec = makeTransparent(col_vec, alpha = 0.35)  
    }
    
    pb = barplot(rbind(colMeans(B_chain[, indices_i] == 1),
                       colMeans(B_chain[, indices_i] == 2),
                       colMeans(B_chain[, indices_i] == 3),
                       colMeans(B_chain[, indices_i] == 4),
                       colMeans(B_chain[, indices_i] == 5)), 
                 col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                 xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
    
    # Heart Rate and MAP double plot -----------------------------------------
    if(mean(use_data[indices_i, 'clinic_rule']) != 0) {
        if(i %in% eid_poor) {
            title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                                mean(use_data[indices_i, 'RBC_rule']),
                                ', clinic = ', mean(use_data[indices_i, 'clinic_rule']),
                                ' (flag)')
        } else {
            title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                                mean(use_data[indices_i, 'RBC_rule']),
                                ', clinic = ', mean(use_data[indices_i, 'clinic_rule']))
        }
    } else {
        if(i %in% eid_poor) {
            title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                                mean(use_data[indices_i, 'RBC_rule']), ' (flag)')   
        } else {
            title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                                mean(use_data[indices_i, 'RBC_rule']))
        }
    }
    
    hr_upper = colQuantiles( Hr_chain[, indices_i, drop=F], probs=.975)
    hr_lower = colQuantiles( Hr_chain[, indices_i, drop=F], probs=.025)
    bp_upper = colQuantiles( Map_chain[, indices_i, drop=F], probs=.975)
    bp_lower = colQuantiles( Map_chain[, indices_i, drop=F], probs=.025)
    
    hr_map_ylim = c(min(hr_lower, bp_lower), max(hr_upper, bp_upper))
    
    # Make a new plot to add the background color
    plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
         xlab='time', ylab=NA, xaxt='n', col.main='green',
         col.axis='green')
    
    rect(xleft = rect_coords$t[-nrow(rect_coords)], 
         ybottom = hr_map_ylim[1], 
         xright = rect_coords$t[-1], 
         ytop = hr_map_ylim[2],
         col = col_vec[-nrow(rect_coords)],
         border = NA)
    
    plotCI( x = pb, y=colMeans(Hr_chain[, indices_i, drop=F]), ui=hr_upper, li=hr_lower,
            main=title_name,
            xlab='time', ylab=NA, xaxt='n', col.main='green',
            col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
            xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim, add =T) 
    plotCI( x = pb, y=colMeans(Map_chain[, indices_i, drop=F]), ui=bp_upper, li=bp_lower,
            main=title_name,
            xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
            col = 'orange',
            xlim = range(pb) + c(-0.5,0.5), add = T) 
    legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'HR', 'MAP'), pch=15, pt.cex=1.5, 
            col=c( 'aquamarine', 'orange'))
    grid( nx=20, NULL, col='white')
    axis( side=1, at=pb, col.axis='green', labels=t_grid)
    
    abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
    abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    
    
    # Hemoglobin and Lactate double plot -------------------------------------
    if(mean(use_data[indices_i, 'clinic_rule']) != 0) {
        title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ', 
                            mean(use_data[indices_i, 'RBC_rule']),
                            ', clinic = ', mean(use_data[indices_i, 'clinic_rule']))
    } else {
        title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ',
                            mean(use_data[indices_i, 'RBC_rule']))
    }
    
    hc_upper = colQuantiles( Hc_chain[, indices_i, drop=F], probs=.975)
    hc_lower = colQuantiles( Hc_chain[, indices_i, drop=F], probs=.025)
    la_upper = colQuantiles( La_chain[, indices_i, drop=F], probs=.975)
    la_lower = colQuantiles( La_chain[, indices_i, drop=F], probs=.025)
    
    hr_map_ylim = c(min(hc_lower, la_lower), max(hc_upper, la_upper))
    
    plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
         xlab='time', ylab=NA, xaxt='n', col.main='green',
         col.axis='green')
    
    rect(xleft = rect_coords$t[-nrow(rect_coords)], 
         ybottom = hr_map_ylim[1], 
         xright = rect_coords$t[-1], 
         ytop = hr_map_ylim[2],
         col = col_vec[-nrow(rect_coords)],
         border = NA)
    
    plotCI(x = pb, y = colMeans(Hc_chain[, indices_i, drop=F]), ui=hc_upper, li=hc_lower,
           main=title_name,
           xlab='time', ylab=NA, xaxt='n', col.main='green',
           col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
           xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim, add = T) 
    plotCI( x = pb, y=colMeans(La_chain[, indices_i, drop=F]), ui=la_upper, li=la_lower,
            main=title_name,
            xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
            col = 'orange',
            xlim = range(pb) + c(-0.5,0.5), add = T) 
    legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'hemo', 'lactate'), pch=15, pt.cex=1.5, 
            col=c( 'aquamarine', 'orange'))
    grid( nx=20, NULL, col='white')
    axis( side=1, at=pb, col.axis='green', labels=t_grid)
    
    abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
    abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    
    # Medication admin plot ----------------------------------------------------
    med_i = Dn_omega[[which(EIDs == i)]]
    omega_i = omega_i_mat[[which(EIDs == i)]]
    med_i_mat = stacked_chains = do.call( rbind, med_i)
    
    hr_med_i_mat = med_i_mat[seq(2, nrow(med_i_mat), by = 4), ]
    map_med_i_mat = med_i_mat[seq(3, nrow(med_i_mat), by = 4), ]
    
    hr_mean_effect = hr_med_i_mat %*% omega_i
    map_mean_effect = map_med_i_mat %*% omega_i
    
    hr_med_i_mat = hr_med_i_mat[, hr_map_names %in% c('hr_cont', 'hr_disc')]
    map_med_i_mat = map_med_i_mat[, hr_map_names %in% c('map_cont', 'map_disc')]
    
    upp_down = c(-1, -1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1,  1,
                 -1,  1,  1,  1, -1,  1, -1,  1, -1, -1, -1, -1, -1,
                 -1, -1, -1, -1,  1, -1,  1, -1, -1,  1,  1, -1,  1,
                 -1, -1,  1,  1, -1, -1, -1, -1, -1, -1,  1, -1,  1,
                 1, -1,  1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,
                 -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, -1,
                 -1, -1,  1, -1, -1, -1, -1, -1, -1,  1)
    
    hr_upp_down = upp_down[hr_map_names %in% c('hr_cont', 'hr_disc')]
    map_upp_down = upp_down[hr_map_names %in% c('map_cont', 'map_disc')]
    
    hr_upp_i   = hr_med_i_mat[,hr_upp_down == 1]
    hr_down_i  = hr_med_i_mat[,hr_upp_down == -1]
    map_upp_i  = map_med_i_mat[,map_upp_down == 1]
    map_down_i = map_med_i_mat[,map_upp_down == -1]
    
    total_hr_up  = rowSums(hr_upp_i)
    total_hr_dn  = rowSums(hr_down_i)
    total_map_up = rowSums(map_upp_i)
    total_map_dn = rowSums(map_down_i)
    
    hr_map_ylim = c(min(total_hr_up, total_hr_dn, total_map_up, total_map_dn,
                        hr_mean_effect, map_mean_effect), 
                    max(total_hr_up, total_hr_dn, total_map_up, total_map_dn,
                        hr_mean_effect, map_mean_effect))
    if(hr_map_ylim[1] == hr_map_ylim[2]) hr_map_ylim = c(0,1)
    
    plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main='Med. admin',
         xlab='time', ylab=NA, xaxt='n', col.main='green',
         col.axis='green')
    
    if(simulation) {
        rect(xleft = rect_coords$t[-nrow(rect_coords)], 
             ybottom = hr_map_ylim[1], 
             xright = rect_coords$t[-1], 
             ytop = hr_map_ylim[2],
             col = col_vec[-nrow(rect_coords)],
             border = NA)
    } 
    points(x = pb, y = hr_mean_effect, xlab='time', ylab=NA, 
           col.main='green', col.axis='green', 
           col = 'aquamarine', pch = 16) 
    points(x = pb, y = map_mean_effect, xlab='time', ylab=NA, 
           col = 'orange', pch = 16) 
    lines(x = pb, y = hr_mean_effect, xlab='time', ylab=NA, 
          lwd=2, lty = 1, col = 'aquamarine') 
    lines(x = pb, y = map_mean_effect, xlab='time', ylab=NA, 
          lwd=2, lty = 1, col = 'orange') 
    
    lines(x = pb, y = total_hr_up, xlab='time', ylab=NA, 
          lwd=1, lty = 2, col = 'aquamarine4') 
    lines(x = pb, y = total_map_up, xlab='time', ylab=NA,
          lwd=1, lty = 3, col = 'darkolivegreen2') 
    lines(x = pb, y = total_hr_dn, xlab='time', ylab=NA,
          lwd=1, lty = 4, col = 'deeppink')
    lines(x = pb, y = total_map_dn, xlab='time', ylab=NA,
          lwd=1, lty = 5, col = 'palevioletred')
    legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'HR effect', 'MAP effect'), pch=15, pt.cex=1.5, 
            col=c( 'aquamarine', 'orange'))
    axis( side=1, at=pb, col.axis='green', labels=t_grid)
    
    abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
    abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    
    # BAR PLOTS --------------------------------------------------------------
    barplot(rbind(colMeans(B_chain[, indices_i] == 1),
                  colMeans(B_chain[, indices_i] == 2),
                  colMeans(B_chain[, indices_i] == 3),
                  colMeans(B_chain[, indices_i] == 4),
                  colMeans(B_chain[, indices_i] == 5)), 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
            xlab='time', space=0, col.main='green', border=NA,
            xlim=range(pb) + c(-0.5,0.5)) 
    grid( nx=20, NULL, col='white')
    legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'Baseline', 'State 2', 'State 3', 'State 4', 'State 5'), 
            pch=15, pt.cex=1.5, 
            col=c( 'dodgerblue', 'firebrick1', 'yellow2','green', 'darkgray'))
    legend( 'topleft', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'RBC order', 'RBC admin'), pch=15, pt.cex=1.5, 
            col=c( 'darkorchid1', 'aquamarine'))				
    axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
    axis( side=2, at=0:1, col.axis='green')
    
    
    abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
    abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
}
dev.off()