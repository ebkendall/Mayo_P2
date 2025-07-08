args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])
simulation = as.logical(as.numeric(args[2]))

if(simulation) {
    it_seq = 1:1
    index_seeds = c(1:50)
} else {
    it_seq = 6:10
    index_seeds = c(1:5)
}

trialNum = 1
states_per_step = 0
steps_per_it = 1

load('Data_sim/true_pars.rda')
load('Data_real/hr_map_names.rda')
load('Data_real/Dn_omega_names.rda')

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

# Loading parameter info -------------------------------------------------------
a_w_chain_id = c(1, 86, 163, 237, 427)
med_check_inds = c(2, 6, 9, 12, 29, 42, 46, 74)

chain_list = vector(mode = "list", length = length(index_seeds))
a_chain_list = vector(mode = 'list', length = length(a_w_chain_id))
w_chain_list = vector(mode = 'list', length = length(a_w_chain_id))

for(a in 1:length(a_chain_list)) {
    a_chain_list[[a]] = vector(mode = 'list', length = length(index_seeds))
    w_chain_list[[a]] = vector(mode = 'list', length = length(index_seeds))
}

ind = 0
post_means_mat = NULL
covg = NULL
for(seed in index_seeds){
    for(it in it_seq) {
        file_name = NULL
        if(simulation) {
            file_name = paste0('Model_out/mcmc_out_', trialNum, '_', seed, 'it', 
                               it, '_samp', sampling_num, '_sim.rda')
        } else {
            file_name = paste0('Model_out/mcmc_out_', trialNum, '_', seed, 'it', 
                               it, '_samp', sampling_num, '_', states_per_step,
                               '_', steps_per_it,'.rda')
        }
        
        if(file.exists(file_name)) {
            load(file_name)
            
            par_index = mcmc_out_temp$par_index
            
            if(which(it_seq == it) == 1) {
                ind = ind + 1
                
                chain_list[[ind]] = mcmc_out_temp$chain
                
                for(a in 1:length(a_chain_list)) {
                    a_chain_list[[a]][[ind]] = mcmc_out_temp$A_chain[[a]]
                    w_chain_list[[a]][[ind]] = mcmc_out_temp$W_chain[[a]]
                }
            } else {
                chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out_temp$chain)
                
                for(a in 1:length(a_chain_list)) {
                    a_chain_list[[a]][[ind]] = cbind(a_chain_list[[a]][[ind]], mcmc_out_temp$A_chain[[a]])
                    w_chain_list[[a]][[ind]] = cbind(w_chain_list[[a]][[ind]], mcmc_out_temp$W_chain[[a]])
                }
            }
            
            print(paste0(ind, ": ", file_name))
            print("accept")
            print(mcmc_out_temp$accept)
            print("pscale")
            print(mcmc_out_temp$pscale)
            
            rm(mcmc_out_temp)    
        } else {
            print(paste0("Missing! ", file_name))
        }
    }
    
    post_means_temp = matrix(colMeans(chain_list[[ind]]), nrow = 1)
    post_means_mat = rbind(post_means_mat, post_means_temp)
    
    if(simulation) {
        covg_low = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.025)})
        covg_high = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.975)})
        covg_temp = as.numeric(true_pars <= covg_high & true_pars >= covg_low)
        covg = rbind(covg, covg_temp)    
    }
}

if(simulation) {mean_covg = colMeans(covg)}

stacked_chains = do.call( rbind, chain_list)
a_stacked_chains = vector(mode = 'list', length = length(a_w_chain_id))
w_stacked_chains = vector(mode = 'list', length = length(a_w_chain_id))
for(a in 1:length(a_stacked_chains)) {
    a_stacked_chains[[a]] = do.call('cbind', a_chain_list[[a]])
    w_stacked_chains[[a]] = do.call('cbind', w_chain_list[[a]])
}

if(!simulation) {
    par_means = colMeans(stacked_chains)
    save(par_means, file = paste0('Model_out/post_means_samp', sampling_num, 
                                '_it', max(it_seq), '_real.rda'))
}

pdf_file = NULL
if(simulation) {
    pdf_file = paste0('Plots/sim_trace_', trialNum, '_samp', sampling_num, '_it', max(it_seq), '.pdf')
} else {
    pdf_file = paste0('Plots/real_trace_', trialNum, '_samp', sampling_num, '_it', max(it_seq), '.pdf')
}

pdf(pdf_file)
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
        
        # Trace plots ----------------------------------------------------------
        x_title_1 = NULL
        if(simulation) {
            x_title_1 = paste0("95% CI: [", round(lower, 4), ", ", round(upper, 4), "]",
                               " Covg =", round(mean_covg[r], 3))
        } else {
            x_title_1 = paste0("95% CI: [", round(lower, 4), ", ", round(upper, 4), "]")
        }
        
        y_limit = range(stacked_chains[,r])
        plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,nrow(chain_list[[1]])),
              ylim=y_limit, xlab = x_title_1, col.main = title_color)
        
        for(seed in 1:length(chain_list)) {
            lines( chain_list[[seed]][,r], type='l', col=seed)
        }
        
        # Posterior samples ----------------------------------------------------
        x_title_2 = NULL
        if(simulation) {
            x_title_2 = paste0('Mean =',toString(parMean),
                               ' Median =',toString(parMedian),
                               ' True =', true_pars[r])    
        } else {
            x_title_2 = paste0('Mean =',toString(parMean),
                               ' Median =',toString(parMedian))
        }
        
        hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, 
              main=NA, freq=FALSE, xlab=x_title_2)
        abline( v=upper, col='red', lwd=2, lty=2)
        abline( v=lower, col='purple', lwd=2, lty=2)
        if(simulation){ abline( v=true_pars[r], col='green', lwd=2, lty=2) }
    }   
}

# Plot the sampled alpha_i -----------------------------------------------------
alpha_i_lab = labels[par_index$vec_alpha_tilde]
for(s in 1:length(a_w_chain_id)) {
    
    for(l in 1:length(alpha_i_lab)) {
        parMean = round( mean(a_stacked_chains[[s]][l,]), 4)
        parMedian = round( median(a_stacked_chains[[s]][l,]), 4)
        upper = quantile(a_stacked_chains[[s]][l,], prob=.975)
        lower = quantile(a_stacked_chains[[s]][l,], prob=.025)
        y_limit = range(a_stacked_chains[[s]][l,])
        
        plot( NULL, ylab=NA, main=paste0(a_w_chain_id[s], ": ", alpha_i_lab[l]), 
              xlim=c(1,ncol(a_chain_list[[s]][[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"),
              col.main = 'black')
        
        for(seed in 1:length(a_chain_list[[s]])) {
            lines( a_chain_list[[s]][[seed]][l,], type='l', col=seed)
        }
        
        x_label = paste0('Mean =',toString(parMean),
                         ' Median =',toString(parMedian))
        
        hist(a_stacked_chains[[s]][l,], breaks=sqrt(ncol(a_stacked_chains[[s]])), 
             ylab=NA, main=NA, freq=FALSE, xlab=x_label)    
    }
}

# Plot the sampled omega_i -----------------------------------------------------
omega_i_lab = labels[par_index$omega_tilde][med_check_inds]
for(s in 1:length(a_w_chain_id)) {
    
    for(l in 1:length(omega_i_lab)) {
        parMean = round( mean(w_stacked_chains[[s]][l,]), 4)
        parMedian = round( median(w_stacked_chains[[s]][l,]), 4)
        upper = quantile(w_stacked_chains[[s]][l,], prob=.975)
        lower = quantile(w_stacked_chains[[s]][l,], prob=.025)
        y_limit = range(w_stacked_chains[[s]][l,])
        
        plot( NULL, ylab=NA, main=paste0(a_w_chain_id[s], ": ", omega_i_lab[l]),
              xlim=c(1,ncol(w_chain_list[[s]][[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"),
              col.main = 'black')
        
        for(seed in 1:length(w_chain_list[[s]])) {
            lines( w_chain_list[[s]][[seed]][l,], type='l', col=seed)
        }
        
        x_label = paste0('Mean =',toString(parMean),
                         ' Median =',toString(parMedian))
        
        hist(w_stacked_chains[[s]][l,], breaks=sqrt(ncol(w_stacked_chains[[s]])),
             ylab=NA, main=NA, freq=FALSE, xlab=x_label)
    }
}

if(simulation) {    
    par(mfrow=c(3, 3))
    lab_ind = 0
    for(s in names(par_index)){
        
        temp_par = par_index[[s]]
        if(s == "vec_sigma_upsilon") {
            temp_par = temp_par[upsilon_ind]
        }
        
        for(r in temp_par){
            lab_ind = r
            
            boxplot(post_means_mat[,r], main = labels[r],
                    xlab = paste0('95% Covg = ', round(mean(covg[,r]), 4)),
                    ylab = paste0('truth = ', true_pars[r]))
            abline(h = true_pars[r], col = 'red')
        }
    }
}


dev.off()
