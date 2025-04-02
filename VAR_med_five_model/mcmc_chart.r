library(matrixStats)
library(plotrix)

trialNum = 1
args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])

it_num = 1
it_seq = 6:30

states_per_step = 0
steps_per_it = 1
S = 5

seed_list = 1:3

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
hemo_results = vector(mode = 'list', length = length(seed_list))
hr_results = vector(mode = 'list', length = length(seed_list))
map_results = vector(mode = 'list', length = length(seed_list))
lact_results = vector(mode = 'list', length = length(seed_list))

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
                           '_', steps_per_it,'.rda')
        load(file_name)
        print(file_name)
        
        # if(it == 1) {
        #     B_chain   = mcmc_out_temp$B_chain[1000:2000, ]
        #     # Hr_chain  = mcmc_out_temp$hr_chain[1000:2000, ]
        #     # Map_chain = mcmc_out_temp$bp_chain[1000:2000, ]
        #     # Hc_chain  = mcmc_out_temp$hc_chain[1000:2000, ]
        #     # La_chain  = mcmc_out_temp$la_chain[1000:2000, ]
        # } else {
        B_chain   = rbind(B_chain, mcmc_out_temp$B_chain[!is.na(mcmc_out_temp$B_chain[,1]),])
        Hr_chain  = rbind(Hr_chain, mcmc_out_temp$hr_chain[!is.na(mcmc_out_temp$hr_chain[,1]),])
        Map_chain = rbind(Map_chain, mcmc_out_temp$bp_chain[!is.na(mcmc_out_temp$bp_chain[,1]),])
        Hc_chain  = rbind(Hc_chain, mcmc_out_temp$hc_chain[!is.na(mcmc_out_temp$hc_chain[,1]),])
        La_chain  = rbind(La_chain, mcmc_out_temp$la_chain[!is.na(mcmc_out_temp$la_chain[,1]),])
        # }
        rm(mcmc_out_temp)
    }    
    
    state_results[[s]] = matrix(nrow = S+1, ncol = ncol(B_chain))
    for(jj in 1:S) {
        state_results[[s]][jj, ] = apply(B_chain, 2, function(x,jj){sum(x == jj)}, jj)
    }
    state_results[[s]][S+1, ] = apply(B_chain, 2, Mode)
    
    hemo_results[[s]] = Hc_chain
    hr_results[[s]] = Hr_chain
    map_results[[s]] = Map_chain
    lact_results[[s]] = La_chain
    
    rm(B_chain)
    rm(Hc_chain)
    rm(Hr_chain)
    rm(Map_chain)
    rm(La_chain)
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
load('Data_real/Dn_omega.rda')
load('Data_real/data_format_train.rda')
load('Data_sim/hr_map_names.rda')
EIDs = unique(data_format[,'EID'])

# Model evaluation plots -------------------------------------------------------
for(one_chart in seed_list) {
    state_counts_col_sum = colSums(state_results[[one_chart]][1:S,])
    state_proportions = matrix(nrow = S, ncol = ncol(state_results[[one_chart]]))
    for(s in 1:S) {
        state_proportions[s,] = state_results[[one_chart]][s, ] / state_counts_col_sum
    }
    
    pdf(paste0('Plots/chart_', trialNum, '_', one_chart, '_samp', sampling_num, '_',
               states_per_step, '_', steps_per_it, '.pdf'))
    panel_dim = c(4,1)
    inset_dim = c(0,-.18)
    par(mfrow=panel_dim, mar=c(2,4,2,4), bg='black', fg='green')
    for(i in EIDs){
        if(which(EIDs == i) %% 100 == 0) print(which(EIDs == i))
        
        indices_i = (data_format[,'EID']==i)
        n_i = sum(indices_i)
        
        t_grid = round(data_format[indices_i, 'time'] / 60, digits = 3)
        t_grid_bar = 1:length(t_grid)
        rbc_times_bar = which(data_format[data_format[,'EID']==i, 'RBC_ordered'] != 0)
        
        rbc_admin_times_bar = which(data_format[data_format[,'EID']==i, 'RBC_admin'] != 0)   
        
        rbc_times = t_grid[rbc_times_bar]
        rbc_admin_times = t_grid[rbc_admin_times_bar]
        
        pb = barplot(state_proportions[,indices_i], 
                     col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                     xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
        
        # Heart Rate and MAP double plot -----------------------------------------
        if(mean(data_format[indices_i, 'clinic_rule']) != 0) {
            title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                                mean(data_format[indices_i, 'RBC_rule']),
                                ', clinic = ', mean(data_format[indices_i, 'clinic_rule']))
        } else {
            title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                                mean(data_format[indices_i, 'RBC_rule']))
        }
        
        hr_upper = colQuantiles( hr_results[[one_chart]][, indices_i, drop=F], probs=.975)
        hr_lower = colQuantiles( hr_results[[one_chart]][, indices_i, drop=F], probs=.025)
        bp_upper = colQuantiles( map_results[[one_chart]][, indices_i, drop=F], probs=.975)
        bp_lower = colQuantiles( map_results[[one_chart]][, indices_i, drop=F], probs=.025)
        
        hr_map_ylim = c(min(hr_lower, bp_lower), max(hr_upper, bp_upper))
        
        # Make a new plot to add the background color
        plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
             xlab='time', ylab=NA, xaxt='n', col.main='green',
             col.axis='green')
        
        plotCI( x = pb, y=colMeans(hr_results[[one_chart]][, indices_i, drop=F]), ui=hr_upper, li=hr_lower,
                main=title_name,
                xlab='time', ylab=NA, xaxt='n', col.main='green',
                col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
                xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim, add =T) 
        plotCI( x = pb, y=colMeans(map_results[[one_chart]][, indices_i, drop=F]), ui=bp_upper, li=bp_lower,
                main=title_name,
                xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
                col = 'orange',
                xlim = range(pb) + c(-0.5,0.5), add = T) 
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'HR', 'MAP'), pch=15, pt.cex=1.5, 
                col=c( 'aquamarine', 'orange'))
        grid( nx=20, NULL, col='white')
        axis( side=1, at=pb, col.axis='green', labels=t_grid)
        
        # abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        
        
        # Hemoglobin and Lactate double plot -------------------------------------
        if(mean(data_format[indices_i, 'clinic_rule']) != 0) {
            title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ', 
                                mean(data_format[indices_i, 'RBC_rule']),
                                ', clinic = ', mean(data_format[indices_i, 'clinic_rule']))
        } else {
            title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ',
                                mean(data_format[indices_i, 'RBC_rule']))
        }
        
        hc_upper = colQuantiles( hemo_results[[one_chart]][, indices_i, drop=F], probs=.975)
        hc_lower = colQuantiles( hemo_results[[one_chart]][, indices_i, drop=F], probs=.025)
        la_upper = colQuantiles( lact_results[[one_chart]][, indices_i, drop=F], probs=.975)
        la_lower = colQuantiles( lact_results[[one_chart]][, indices_i, drop=F], probs=.025)
        
        hr_map_ylim = c(min(hc_lower, la_lower), max(hc_upper, la_upper))
        
        plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
             xlab='time', ylab=NA, xaxt='n', col.main='green',
             col.axis='green')
        
        plotCI(x = pb, y = colMeans(hemo_results[[one_chart]][, indices_i, drop=F]), ui=hc_upper, li=hc_lower,
               main=title_name,
               xlab='time', ylab=NA, xaxt='n', col.main='green',
               col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
               xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim, add = T) 
        plotCI( x = pb, y=colMeans(lact_results[[one_chart]][, indices_i, drop=F]), ui=la_upper, li=la_lower,
                main=title_name,
                xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
                col = 'orange',
                xlim = range(pb) + c(-0.5,0.5), add = T) 
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'hemo', 'lactate'), pch=15, pt.cex=1.5, 
                col=c( 'aquamarine', 'orange'))
        grid( nx=20, NULL, col='white')
        axis( side=1, at=pb, col.axis='green', labels=t_grid)
        
        # abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        
        # # Medication admin plot ----------------------------------------------------
        # med_i = Dn_omega[[which(EIDs == i)]]
        # omega_i = omega_i_mat[[which(EIDs == i)]]
        # med_i_mat = stacked_chains = do.call( rbind, med_i)
        # 
        # hr_med_i_mat = med_i_mat[seq(2, nrow(med_i_mat), by = 4), ]
        # map_med_i_mat = med_i_mat[seq(3, nrow(med_i_mat), by = 4), ]
        # 
        # hr_mean_effect = hr_med_i_mat %*% omega_i
        # map_mean_effect = map_med_i_mat %*% omega_i
        # 
        # hr_med_i_mat = hr_med_i_mat[, hr_map_names %in% c('hr_cont', 'hr_disc')]
        # map_med_i_mat = map_med_i_mat[, hr_map_names %in% c('map_cont', 'map_disc')]
        # 
        # upp_down = c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
        #              -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
        #              -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
        #              -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
        # 
        # hr_upp_down = upp_down[hr_map_names %in% c('hr_cont', 'hr_disc')]
        # map_upp_down = upp_down[hr_map_names %in% c('map_cont', 'map_disc')]
        # 
        # hr_upp_i   = hr_med_i_mat[,hr_upp_down == 1]
        # hr_down_i  = hr_med_i_mat[,hr_upp_down == -1]
        # map_upp_i  = map_med_i_mat[,map_upp_down == 1]
        # map_down_i = map_med_i_mat[,map_upp_down == -1]
        # 
        # total_hr_up  = rowSums(hr_upp_i)
        # total_hr_dn  = rowSums(hr_down_i)
        # total_map_up = rowSums(map_upp_i)
        # total_map_dn = rowSums(map_down_i)
        # 
        # hr_map_ylim = c(min(total_hr_up, total_hr_dn, total_map_up, total_map_dn,
        #                     hr_mean_effect, map_mean_effect), 
        #                 max(total_hr_up, total_hr_dn, total_map_up, total_map_dn,
        #                     hr_mean_effect, map_mean_effect))
        # if(hr_map_ylim[1] == hr_map_ylim[2]) hr_map_ylim = c(0,1)
        
        plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main='Med. admin',
             xlab='time', ylab=NA, xaxt='n', col.main='green',
             col.axis='green')
        
        # points(x = pb, y = hr_mean_effect, xlab='time', ylab=NA, 
        #        col.main='green', col.axis='green', 
        #        col = 'aquamarine', pch = 16) 
        # points(x = pb, y = map_mean_effect, xlab='time', ylab=NA, 
        #        col = 'orange', pch = 16) 
        # lines(x = pb, y = hr_mean_effect, xlab='time', ylab=NA, 
        #       lwd=2, lty = 1, col = 'aquamarine') 
        # lines(x = pb, y = map_mean_effect, xlab='time', ylab=NA, 
        #       lwd=2, lty = 1, col = 'orange') 
        # 
        # lines(x = pb, y = total_hr_up, xlab='time', ylab=NA, 
        #       lwd=1, lty = 2, col = 'aquamarine4') 
        # lines(x = pb, y = total_map_up, xlab='time', ylab=NA,
        #       lwd=1, lty = 3, col = 'darkolivegreen2') 
        # lines(x = pb, y = total_hr_dn, xlab='time', ylab=NA,
        #       lwd=1, lty = 4, col = 'deeppink')
        # lines(x = pb, y = total_map_dn, xlab='time', ylab=NA,
        #       lwd=1, lty = 5, col = 'palevioletred')
        # legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
        #         legend=c( 'HR effect', 'MAP effect'), pch=15, pt.cex=1.5, 
        #         col=c( 'aquamarine', 'orange'))
        # axis( side=1, at=pb, col.axis='green', labels=t_grid)
        # 
        # # abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        # abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        
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
        
        
        # abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    }
    dev.off()
}

