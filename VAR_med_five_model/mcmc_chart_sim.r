library(matrixStats)
library(plotrix)

trialNum = 1
args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])
one_chart = as.numeric(args[2])

it_num = 1
it_seq = 1:it_num

states_per_step = 0
steps_per_it = 1
S = 5

seed_list = c(1:5)

df_num = 1

# Mode of the state sequences -------------------------------------------------
Mode <- function(x) {
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
}

# Function to change transparency of colors # ----------------------------------
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

# Load the model output -------------------------------------------------------
state_results = vector(mode = 'list', length = length(seed_list))

for(s in 1:length(seed_list)) {
    
    B_chain   = NULL
    Hr_chain  = NULL
    Map_chain = NULL
    Hc_chain  = NULL
    La_chain  = NULL
    
    seed_num = seed_list[s]
    for(it in it_seq) {
        
        file_name = paste0('Model_out/mcmc_out_', trialNum,'_', seed_num,
                           'it', it, '_samp', sampling_num, '_', states_per_step,
                           '_', steps_per_it,'_sim.rda')
        load(file_name)
        print(file_name)
        
        keep_ind = seq(1, nrow(mcmc_out_temp$B_chain), by = 5)
        # if(it == 1) {
        #     B_chain   = mcmc_out_temp$B_chain[1000:2000, ]
        #     # Hr_chain  = mcmc_out_temp$hr_chain[1000:2000, ]
        #     # Map_chain = mcmc_out_temp$bp_chain[1000:2000, ]
        #     # Hc_chain  = mcmc_out_temp$hc_chain[1000:2000, ]
        #     # La_chain  = mcmc_out_temp$la_chain[1000:2000, ]
        # } else {
            B_chain   = rbind(B_chain, mcmc_out_temp$B_chain[keep_ind, ])
            # Hr_chain  = rbind(Hr_chain, mcmc_out_temp$hr_chain)
            # Map_chain = rbind(Map_chain, mcmc_out_temp$bp_chain)
            # Hc_chain  = rbind(Hc_chain, mcmc_out_temp$hc_chain)
            # La_chain  = rbind(La_chain, mcmc_out_temp$la_chain)
        # }
        rm(mcmc_out_temp)
    }    
    
    state_results[[s]] = matrix(nrow = S+1, ncol = ncol(B_chain))
    for(jj in 1:S) {
        state_results[[s]][jj, ] = apply(B_chain, 2, function(x,jj){sum(x == jj)}, jj)
    }
    state_results[[s]][S+1, ] = apply(B_chain, 2, Mode)
    rm(B_chain)
}

# Summarize data into combo and indiv results ----------------------------------
combo_counts = state_results[[1]]
for(i in 2:length(state_results)) {
    combo_counts = combo_counts + state_results[[i]]
}
combo_counts[S+1, ] = apply(combo_counts[1:S,], 2, which.max)

combo_counts_col_sum = colSums(combo_counts[1:S,])
combo_counts_props = matrix(nrow = S, ncol = ncol(combo_counts))
for(s in 1:S) {
    combo_counts_props[s,] = combo_counts[s, ] / combo_counts_col_sum
}

all_seeds_state_mode = matrix(nrow = length(seed_list)+1, ncol = ncol(combo_counts))
for(i in 1:length(seed_list)) {
    all_seeds_state_mode[i, ] = state_results[[i]][S+1,]
}
all_seeds_state_mode[length(seed_list)+1, ] = apply(combo_counts, 2, which.max)

# Load simulated data ----------------------------------------------------------
load(paste0('Data_sim/use_data_', df_num, '.rda'))
load(paste0('Data_sim/alpha_i_mat_', df_num, '.rda'))
load(paste0('Data_sim/omega_i_mat_', df_num, '.rda'))
load('Data_sim/Dn_omega_sim.rda')
Dn_omega = Dn_omega_sim
rm(Dn_omega_sim)
EIDs = unique(use_data[,'EID'])

ss_truth = use_data[,"b_true"]

Hr_chain  = t(use_data[,"hr",drop=F])
Map_chain = t(use_data[,"map",drop=F])
Hc_chain  = t(use_data[,"hemo",drop=F])
La_chain  = t(use_data[,"lactate",drop=F])

# Summary stats of state identification ----------------------------------------
print("Summary of identifying correct states with mode")
for(s in 1:(length(seed_list)+1)) {
    if(s <= length(seed_list)) {
        print(paste0("Seed ", s))
    } else {
        print("Combo")
    }
    
    if(length(ss_truth) != length(all_seeds_state_mode[s, ])) {
        print("ERROR")
    } else {
        print(sum(ss_truth == all_seeds_state_mode[s, ]) / length(all_seeds_state_mode[s, ]))   
    }
}

print("Average (across subject) of proportion of correct states with mode")
mode_correctness = matrix(nrow = length(EIDs), ncol = length(seed_list)+1)
sensitivity_S2 = matrix(nrow = length(EIDs), ncol = length(seed_list)+1)
specificity_S2 = matrix(nrow = length(EIDs), ncol = length(seed_list)+1)
for(s in 1:(length(seed_list)+1)) {
    
    prop_sub = rep(0, length(EIDs))

    for(j in 1:length(EIDs)) {
        sub_ind_j = which(use_data[,"EID"] == EIDs[j])
        ss_truth_j = ss_truth[sub_ind_j]
        state_seq_mode_j = all_seeds_state_mode[s, sub_ind_j]
        
        prop_sub[j] = sum(ss_truth_j == state_seq_mode_j) / length(state_seq_mode_j)
    }
    
    if(s <= length(seed_list)) {
        print(paste0("Seed ", s))
    } else {
        print("Combo")
    }
    
    print(summary(prop_sub))
    mode_correctness[,s] = prop_sub
    
    # eid_poor = NULL
    # eid_poor = EIDs[prop_sub < 0.9]
    
    # Sensitivity of state 2: Pr(predict S2 | true S2)
    predict_at_true_S2 = all_seeds_state_mode[s, (ss_truth == 2)]
    print(paste0("Sensitivity of S2 = ", mean(predict_at_true_S2 == 2)))
    
    # Specificity of state 2: Pr(predict not S2 | true not S2)
    predict_not_S2 =  all_seeds_state_mode[s, (ss_truth != 2)]
    print(paste0("Specificity of S2 = ", mean(predict_not_S2 != 2)))
}

# Choose a subset of the subjects to plot --------------------------------------
EID_plot = unique(c(use_data[use_data[,"RBC_rule"] != 0,"EID"],
                    use_data[use_data[,"clinic_rule"] != 0,"EID"]))
EID_not_chosen_yet = EIDs[!(EIDs %in% EID_plot)]
set.seed(2025)
EID_plot = c(EID_plot, sample(x = EID_not_chosen_yet, 
                              size = 200 - length(EID_plot), 
                              replace = F))

# Model evaluation plots -------------------------------------------------------
state_proportions = NULL
if(one_chart == 0) {
    state_proportions = combo_counts_props
} else {
    state_counts_col_sum = colSums(state_results[[one_chart]][1:S,])
    state_proportions = matrix(nrow = S, ncol = ncol(state_results[[one_chart]]))
    for(s in 1:S) {
        state_proportions[s,] = state_results[[one_chart]][s, ] / state_counts_col_sum
    }
}

load('Data_sim/hr_map_names.rda')

pdf(paste0('Plots/sim_chart_', trialNum, '_', one_chart, '_samp', 
           sampling_num, '_it', max(it_seq), '.pdf'))
panel_dim = c(4,1)
inset_dim = c(0,-.18)
par(mfrow=panel_dim, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EID_plot){
    
    if(which(EID_plot == i) %% 100 == 0) print(which(EID_plot == i))
    
    indices_i = (use_data[,'EID']==i)
    n_i = sum(indices_i)
    
    t_grid = round(use_data[indices_i, 'time'] / 60, digits = 3)
    t_grid_bar = 1:length(t_grid)
    rbc_times_bar = which(use_data[use_data[,'EID']==i, 'RBC_ordered'] != 0)
    
    rbc_admin = c(head(use_data[use_data[,'EID']==i, "n_RBC_admin"], 1),
                  diff(use_data[use_data[,'EID']==i, "n_RBC_admin"]))
    rbc_admin_times_bar = which(rbc_admin != 0)
    
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
    
    pb = barplot(state_proportions[,indices_i], 
                 col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                 xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
    
    # Heart Rate and MAP double plot -----------------------------------------
    if(mean(use_data[indices_i, 'clinic_rule']) != 0) {
        title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                            mean(use_data[indices_i, 'RBC_rule']),
                            ', clinic = ', mean(use_data[indices_i, 'clinic_rule']))
    } else {
        title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                            mean(use_data[indices_i, 'RBC_rule']))
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
    med_i_mat = do.call( rbind, med_i)
    
    hr_med_i_mat = med_i_mat[seq(2, nrow(med_i_mat), by = 4), ]
    map_med_i_mat = med_i_mat[seq(3, nrow(med_i_mat), by = 4), ]
    
    hr_mean_effect = hr_med_i_mat %*% omega_i
    map_mean_effect = map_med_i_mat %*% omega_i
    
    hr_med_i_mat = hr_med_i_mat[, hr_map_names %in% c('hr_cont', 'hr_disc')]
    map_med_i_mat = map_med_i_mat[, hr_map_names %in% c('map_cont', 'map_disc')]
    
    upp_down = c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                 -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                 -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                 -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
    
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
    
    rect(xleft = rect_coords$t[-nrow(rect_coords)], 
         ybottom = hr_map_ylim[1], 
         xright = rect_coords$t[-1], 
         ytop = hr_map_ylim[2],
         col = col_vec[-nrow(rect_coords)],
         border = NA)
    
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
    barplot(state_proportions[,indices_i], 
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
