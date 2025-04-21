args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])

trialNum = 1
it_seq = 1:5

# Loading parameter info -------------------------------------------------------
load(paste0('Model_out/mcmc_results_', trialNum, '_samp', 
            sampling_num, '_it', max(it_seq), '.rda'))
load(paste0('Model_out/re_alpha_', trialNum, '_samp', 
            sampling_num, '_it', max(it_seq), '.rda'))
# load(paste0('Model_out/re_omega_', trialNum, '_samp', 
#             sampling_num, '_it', max(it_seq), '.rda'))

load('Data_sim/hr_map_names.rda')
load('Data_real/Dn_omega_names.rda')

stacked_chains = mcmc_results[[1]]$par
for(s in 2:length(mcmc_results)) {
    stacked_chains = rbind(stacked_chains, mcmc_results[[s]]$par)
}

# Loading random effect --------------------------------------------------------
a_w_chain_id = c(1, 86, 163, 237, 427)
med_check_inds = c(2, 6, 9, 12, 29, 42, 46, 74)

a_stacked_chains = vector(mode = 'list', length = length(rand_eff_alpha))
# w_stacked_chains = vector(mode = 'list', length = length(rand_eff_omega))
for(a in 1:length(a_stacked_chains)) {
    a_stacked_chains[[a]] = do.call('cbind', rand_eff_alpha[[a]])
    # w_stacked_chains[[a]] = do.call('cbind', rand_eff_omega[[a]])
}


# Labels -----------------------------------------------------------------------
labels = c("beta (n_RBC_admin): hemo", "beta (n_RBC_admin): hr", 
           "beta (n_RBC_admin): map", "beta (n_RBC_admin): lact",
           "intercept (hemo)", "slope bleeding (hemo)", 
           "slope recovery (hemo)", "slope state 4 (hemo)", "slope state 5 (hemo)",
           "intercept (hr)", "slope bleeding (hr)", 
           "slope recovery (hr)", "slope state 4 (hr)", "slope state 5 (hr)",
           "intercept (map)", "slope bleeding (map)", 
           "slope recovery (map)", "slope state 4 (map)", "slope state 5 (map)",
           "intercept (lact)", "slope bleeding (lact)", 
           "slope recovery (lact)", "slope state 4 (lact)", "slope state 5 (lact)",
           paste0("Upsilon (", 1:20, ", ", rep(1:20, each = 20), ")"), 
           "AR (hemoglobin)", "AR (heart rate)", "AR (MAP)", "AR (lactate)",
           "Var(hemo)", "Cov(hemo, hr)", "Cov(hemo, map)", "Cov(hemo, lact)", 
           "Cov(hr, hemo)", "Var(hr)", "Cov(hr, map)", "Cov(hr, lact)",
           "Cov(map, hemo)", "Cov(map, hr)", "Var(map)", "Cov(map, lact)",
           "Cov(lact, hemo)", "Cov(lact, hr)", "Cov(lact, map)", "Var(lact)",
           "intercept: S1 --> S2", "RBC_order: S1 --> S2",  
           "intercept: S1 --> S4", "RBC_order: S1 --> S4",
           "intercept: S2 --> S3", "RBC_order: S2 --> S3",  
           "intercept: S2 --> S4", "RBC_order: S2 --> S4", 
           "intercept: S3 --> S1", "RBC_order: S3 --> S1",  
           "intercept: S3 --> S2", "RBC_order: S3 --> S2",
           "intercept: S3 --> S4", "RBC_order: S3 --> S4",  
           "intercept: S4 --> S2", "RBC_order: S4 --> S2",
           "intercept: S4 --> S5", "RBC_order: S4 --> S5",  
           "intercept: S5 --> S1", "RBC_order: S5 --> S1",
           "intercept: S5 --> S2", "RBC_order: S5 --> S2",  
           "intercept: S5 --> S4", "RBC_order: S5 --> S4",
           "logit Pr(init S2)", "logit Pr(init S3)",
           "logit Pr(init S4)", "logit Pr(init S5)",
           paste0("mean (hr): ", Dn_omega_names[1:34]), 
           paste0("mean (map): ", Dn_omega_names[35:84]), 
           paste0("log Upsilon (hr): ", Dn_omega_names[1:34]), 
           paste0("log Upsilon (map): ", Dn_omega_names[35:84])) 

# Indexing initialization ------------------------------------------------------
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:24
par_index$vec_sigma_upsilon = 25:424
par_index$vec_A = 425:428
par_index$vec_R = 429:444
par_index$vec_zeta = 445:468
par_index$vec_init = 469:472
par_index$omega_tilde = 473:556
par_index$vec_upsilon_omega = 557:640

upsilon_ind = matrix(1:length(par_index$vec_sigma_upsilon), ncol = 20)
upsilon_ind[upper.tri(upsilon_ind, diag = F)] = 0
upsilon_ind = c(upsilon_ind)
upsilon_ind = upsilon_ind[upsilon_ind != 0]

post_means = colMeans(stacked_chains)
save(post_means, file = paste0('Model_out/post_means_samp', sampling_num, 
                               '_it', max(it_seq), '_real.rda'))

# Plotting ---------------------------------------------------------------------
pdf(paste0('Plots/real_trace_', trialNum, '_samp', sampling_num, '_it', max(it_seq), '.pdf'))
par(mfrow=c(3, 2))
lab_ind = 0
for(s in names(par_index)){
    
    temp_par = par_index[[s]]
    if(s == "vec_sigma_upsilon") {
        temp_par = temp_par[upsilon_ind]
    }
    
    for(r in temp_par){
        print(r)
        lab_ind = r
        parMean = round( mean(stacked_chains[,r]), 4)
        parMedian = round( median(stacked_chains[,r]), 4)
        upper = quantile( stacked_chains[,r], prob=.975)
        lower = quantile( stacked_chains[,r], prob=.025)
        
        title_color = "black"
        if(s == "omega_tilde") {
            if(0 < lower) {
                title_color = "red"
            }
            if(0 > upper) {
                title_color = "red"
            }
        }
        
        y_limit = range(stacked_chains[,r])
        plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,nrow(mcmc_results[[1]]$par)),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"),
              col.main = title_color)
        
        for(seed in 1:length(mcmc_results)) {
            lines( mcmc_results[[seed]]$par[,r], type='l', col=seed)
        }
        
        x_label = paste0('Mean =',toString(parMean),
                         ' Median =',toString(parMedian))
        
        hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA, freq=FALSE,
              xlab=x_label)
        abline( v=upper, col='red', lwd=2, lty=2)
        abline( v=lower, col='purple', lwd=2, lty=2)
    }   
}

# Plot the sampled alpha_i -----------------------------------------------------
alpha_i_lab = labels[par_index$vec_alpha_tilde]
for(a in 1:length(rand_eff_alpha)) {
    
    for(l in 1:length(alpha_i_lab)) {
        parMean = round( mean(a_stacked_chains[[a]][l,]), 4)
        parMedian = round( median(a_stacked_chains[[a]][l,]), 4)
        upper = quantile(a_stacked_chains[[a]][l,], prob=.975)
        lower = quantile(a_stacked_chains[[a]][l,], prob=.025)
        y_limit = range(a_stacked_chains[[a]][l,])
        
        plot( NULL, ylab=NA, main=paste0(a_w_chain_id[a], ": ", alpha_i_lab[l]), 
              xlim=c(1,ncol(rand_eff_alpha[[a]][[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"),
              col.main = 'black')
        
        for(seed in 1:length(rand_eff_alpha[[a]])) {
            lines( rand_eff_alpha[[a]][[seed]][l,], type='l', col=seed)
            # abline(h = chain_list[[seed]][1,r], col=seed)
        }
        
        x_label = paste0('Mean =',toString(parMean),
                         ' Median =',toString(parMedian))
        
        hist(a_stacked_chains[[a]][l,], breaks=sqrt(ncol(a_stacked_chains[[a]])), 
             ylab=NA, main=NA, freq=FALSE, xlab=x_label)    
    }
}

# # Plot the sampled omega_i -----------------------------------------------------
# omega_i_lab = labels[par_index$omega_tilde][med_check_inds]
# for(a in 1:length(rand_eff_omega)) {
#     
#     for(l in 1:length(omega_i_lab)) {
#         parMean = round( mean(w_stacked_chains[[a]][l,]), 4)
#         parMedian = round( median(w_stacked_chains[[a]][l,]), 4)
#         upper = quantile(w_stacked_chains[[a]][l,], prob=.975)
#         lower = quantile(w_stacked_chains[[a]][l,], prob=.025)
#         y_limit = range(w_stacked_chains[[a]][l,])
#         
#         plot( NULL, ylab=NA, main=paste0(a_w_chain_id[a], ": ", omega_i_lab[l]),
#               xlim=c(1,ncol(rand_eff_omega[[a]][[1]])),
#               ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
#                                           ", ", round(upper, 4), "]"),
#               col.main = 'black')
#         
#         for(seed in 1:length(rand_eff_omega[[a]])) {
#             lines(rand_eff_omega[[a]][[seed]][l,], type='l', col=seed)
#         }
#         
#         x_label = paste0('Mean =',toString(parMean),
#                          ' Median =',toString(parMedian))
#         
#         hist(w_stacked_chains[[a]][l,], breaks=sqrt(ncol(w_stacked_chains[[a]])),
#              ylab=NA, main=NA, freq=FALSE, xlab=x_label)
#     }
# }

dev.off()