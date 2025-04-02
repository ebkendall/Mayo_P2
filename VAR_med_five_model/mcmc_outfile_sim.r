dir = 'Model_out/'

args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])

trialNum = 1
itNum = 1
index_seeds = c(1:3)
states_per_step = 0
steps_per_it = 1
long_chain = T

load('Data_sim/true_pars.rda')

# Size of posterior sample from mcmc chains
steps = 1001

load('../Data_cleaning/Data/Dn_omega_names.rda')
load('../Data_cleaning/Data/hr_map_names.rda')

labels = c("beta (n_RBC_admin): hemo", "beta (n_RBC_admin): hr", 
           "beta (n_RBC_admin): map", "beta (n_RBC_admin): lact",
           "intercept (hemo)", "slope bleeding (hemo)", "slope recovery (hemo)", "slope state 4 (hemo)", "slope state 5 (hemo)",
           "intercept (hr)", "slope bleeding (hr)", "slope recovery (hr)", "slope state 4 (hr)", "slope state 5 (hr)",
           "intercept (map)", "slope bleeding (map)", "slope recovery (map)", "slope state 4 (map)", "slope state 5 (map)",
           "intercept (lact)", "slope bleeding (lact)", "slope recovery (lact)", "slope state 4 (lact)", "slope state 5 (lact)",
           paste0("Upsilon (", 1:20, ", ", rep(1:20, each = 20), ")"), 
           "AR (hemoglobin)", "AR (heart rate)", "AR (MAP)", "AR (lactate)",
           "Var(hemo)", "Cov(hemo, hr)", "Cov(hemo, map)", "Cov(hemo, lact)", 
           "Cov(hr, hemo)", "Var(hr)", "Cov(hr, map)", "Cov(hr, lact)",
           "Cov(map, hemo)", "Cov(map, hr)", "Var(map)", "Cov(map, lact)",
           "Cov(lact, hemo)", "Cov(lact, hr)", "Cov(lact, map)", "Var(lact)",
           "intercept: S1 --> S2", "RBC_order: S1 --> S2",  "intercept: S1 --> S4", "RBC_order: S1 --> S4",
           "intercept: S2 --> S3", "RBC_order: S2 --> S3",  "intercept: S2 --> S4", "RBC_order: S2 --> S4", 
           "intercept: S3 --> S1", "RBC_order: S3 --> S1",  "intercept: S3 --> S2", "RBC_order: S3 --> S2",
           "intercept: S3 --> S4", "RBC_order: S3 --> S4",  "intercept: S4 --> S2", "RBC_order: S4 --> S2",
           "intercept: S4 --> S5", "RBC_order: S4 --> S5",  "intercept: S5 --> S1", "RBC_order: S5 --> S1",
           "intercept: S5 --> S2", "RBC_order: S5 --> S2",  "intercept: S5 --> S4", "RBC_order: S5 --> S4",
           "logit Pr(init S2)", "logit Pr(init S3)","logit Pr(init S4)", "logit Pr(init S5)",
           paste0("mean (hr): ", Dn_omega_names[1:34]), paste0("mean (map): ", Dn_omega_names[35:84]), 
           paste0("log Upsilon (hr): ", Dn_omega_names[1:34]), paste0("log Upsilon (map): ", Dn_omega_names[35:84])) 
additional_labels = c("Gamma(1,1) stable", "Gamma(2,2) stable", "Gamma(3,3) stable", "Gamma(4,4) stable",
                      "Gamma(1,1) bleed", "Gamma(2,2) bleed", "Gamma(3,3) bleed", "Gamma(4,4) bleed",
                      "Gamma(1,1) recov", "Gamma(2,2) recov", "Gamma(3,3) recov", "Gamma(4,4) recov",
                      "Gamma(1,1) state 4", "Gamma(2,2) state 4", "Gamma(3,3) state 4", "Gamma(4,4) state 4",
                      "Gamma(1,1) state 5", "Gamma(2,2) state 5", "Gamma(3,3) state 5", "Gamma(4,4) state 5")

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))
a_chain_list = vector(mode = 'list', length = length(index_seeds))

for(a in 1:length(a_chain_list)) {
    a_chain_list[[a]] = vector(mode = 'list', length = 10)
}

ind = 0

post_means = NULL
covg = NULL

for(seed in index_seeds){
    
    if(long_chain) {
        it_seq = 1:itNum
    } else {
        it_seq = itNum
    }
    
    for(it in it_seq) {
        
        file_name = paste0(dir,'mcmc_out_', trialNum, '_', seed, 'it', 
                           it, '_samp', sampling_num, '_', states_per_step,
                           '_', steps_per_it,'_sim.rda') 
        
        if(file.exists(file_name)) {
            load(file_name)
            print(paste0(ind, ": ", file_name))
            print("accept")
            print(mcmc_out_temp$accept)
            print("pscale")
            print(mcmc_out_temp$pscale)
            
            par_index = mcmc_out_temp$par_index
            
            if(it == 1) {
                ind = ind + 1
                
                chain_list[[ind]] = mcmc_out_temp$chain
                # for(a in 1:length(a_chain_list[[ind]])) {
                #     a_chain_list[[ind]][[a]] = mcmc_out_temp$A_chain[[a]]
                # }
            } else {
                chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out_temp$chain)
                # for(a in 1:length(a_chain_list[[ind]])) {
                #     a_chain_list[[ind]][[a]] = cbind(a_chain_list[[ind]][[a]], mcmc_out_temp$A_chain[[a]])
                # }
            }
            
            post_means_temp = matrix(colMeans(chain_list[[ind]]), nrow = 1)
            post_means = rbind(post_means, post_means_temp)
            
            covg_low = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.025)})
            covg_high = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.975)})
            covg_temp = as.numeric(true_pars <= covg_high & true_pars >= covg_low)
            covg = rbind(covg, covg_temp)
            
            rm(mcmc_out_temp)    
        } else {
            print("Missing!")
            print(paste0(ind, ": ", file_name))
        }
    }
}

stacked_chains = do.call( rbind, chain_list)

# covg_pars = colMeans(covg)
# boxplot(covg_pars)
# mean(covg_pars)
# sd(covg_pars)

upsilon_ind = matrix(1:length(par_index$vec_sigma_upsilon), ncol = 20)
upsilon_ind[upper.tri(upsilon_ind, diag = F)] = 0
upsilon_ind = c(upsilon_ind)
upsilon_ind = upsilon_ind[upsilon_ind != 0]

# # Re-calculating the Upsilon matrix
# true_gamma = NULL
# 
# gamma_chain = matrix(nrow = nrow(stacked_chains), ncol = 20)
# for(i in 1:nrow(stacked_chains)) {
#     R = matrix(stacked_chains[i, par_index$vec_R], ncol = 4)
#     vec_A1 = stacked_chains[i, par_index$vec_A]
#     scale_A1 = exp(vec_A1) / (1+exp(vec_A1))
# 
#     diag_gamma = c(R[1,1] / (scale_A1[1]^2), R[2,2] / (scale_A1[2]^2),
#                    R[3,3] / (scale_A1[3]^2), R[4,4] / (scale_A1[4]^2),
#                    R[1,1] / (scale_A1[5]^2), R[2,2] / (scale_A1[6]^2),
#                    R[3,3] / (scale_A1[7]^2), R[4,4] / (scale_A1[8]^2),
#                    R[1,1] / (scale_A1[9]^2), R[2,2] / (scale_A1[10]^2),
#                    R[3,3] / (scale_A1[11]^2), R[4,4] / (scale_A1[12]^2),
#                    R[1,1] / (scale_A1[13]^2), R[2,2] / (scale_A1[14]^2),
#                    R[3,3] / (scale_A1[15]^2), R[4,4] / (scale_A1[16]^2),
#                    R[1,1] / (scale_A1[17]^2), R[2,2] / (scale_A1[18]^2),
#                    R[3,3] / (scale_A1[19]^2), R[4,4] / (scale_A1[20]^2))
# 
#     gamma_chain[i, ] = diag_gamma
# }

pdf(paste0('Plots/trace_', trialNum, '_samp', sampling_num, '_',
           states_per_step, '_', steps_per_it, '_1000.pdf'))
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
        plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,nrow(chain_list[[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"),
              col.main = title_color)
        
        for(seed in 1:length(chain_list)) {
            lines( chain_list[[seed]][,r], type='l', col=seed)
            # abline(h = chain_list[[seed]][1,r], col=seed)
        }
        
        x_label = paste0('Mean =',toString(parMean),
                         ' Median =',toString(parMedian),
                         ' True =', round(true_pars[r], 3))
        
        hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA, freq=FALSE,
              xlab=x_label)
        abline( v=upper, col='red', lwd=2, lty=2)
        abline( v=lower, col='purple', lwd=2, lty=2)
        abline( v=true_pars[r], col='green', lwd=2, lty=2)
    }   
}

# par(mfrow=c(3, 3))
# lab_ind = 0
# for(s in names(par_index)){
#     
#     temp_par = par_index[[s]]
#     if(s == "vec_sigma_upsilon") {
#         temp_par = temp_par[upsilon_ind]
#     }
#     
#     for(r in temp_par){
#         lab_ind = r
#         
#         boxplot(post_means[,r], main = labels[r],
#                 xlab = paste0('95% Covg = ', round(mean(covg[,r]), 4)))
#         abline(h = true_pars[r], col = 'red')
#     }
# }


dev.off()
