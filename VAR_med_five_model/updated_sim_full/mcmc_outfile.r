index_seeds = 1:10
trialNum = 1
simulation = F

it_num = 10
max_ind = 14
start_ind = 10

# Parameter initialization -----------------------------------------------------
par_index = list()
par_index$beta = 1:4
par_index$alpha_tilde = 5:20
par_index$upsilon = 21:276
par_index$A = 277:280
par_index$R = 281:296
par_index$zeta = 297:320
par_index$init = 321:324
par_index$omega_tilde = 325:408
par_index$eta_omega = 409:492
par_index$G = 493:508

true_par = rep(0, max(do.call('c', par_index)))
true_par[par_index$beta] = c(0.25, -2, 2, -0.25) # one unit of RBC -> 1 unit increase in hemo in 1 hour
true_par[par_index$alpha_tilde] = c( -5,   5, -2,  2,
                                     10, -10,  2, -2,
                                    -10,  10,  2, -2,
                                      5,  -5, -2,  2)
true_par[par_index$upsilon] = c(diag(c(4, 4, 1, 1, 
                                       4, 4, 1, 1, 
                                       4, 4, 1, 1, 
                                       4, 4, 1, 1)))
true_par[par_index$A] = rep(0, 4)
true_par[par_index$R] = c(diag(c(sqrt(3), sqrt(3), sqrt(3), sqrt(3))))
#    transitions:          1->2,         1->4,         2->3,         2->4, 
#                          3->1,         3->2,         3->4,         4->2, 
#                          4->5,         5->1,         5->2,         5->4
true_par[par_index$zeta] = c(-3.7405, 2.5, -4.2152,   1, -2.6473,-0.5, -2.1475, -0.2, 
                             -3.4459,  -1, -2.9404,   1, -3.2151,   1, -3.1778,  1.5, 
                             -2.0523,   0, -3.4459,-0.2, -3.2404, 2.5, -3.2151,    1)
true_par[par_index$init] = c(0, 0, 0, 0)
true_par[par_index$omega_tilde]= 2 * c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                                       -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                                       -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
true_par[par_index$eta_omega] = rep(1, length(par_index$eta_omega))
true_par[par_index$G] = c(diag(c(sqrt(3), sqrt(3), sqrt(3), sqrt(3))))

load('Data/Dn_omega_names.rda')

labels = c("beta (hemo)", "beta (hr)", "beta (map)", "beta (lact)",
           "slope S2 (hemo)", "slope S3 (hemo)", "slope S4 (hemo)", "slope S5 (hemo)",
           "slope S2 (hr)", "slope S3 (hr)", "slope S4 (hr)", "slope S5 (hr)",
           "slope S2 (map)", "slope S3 (map)", "slope S4 (map)", "slope S5 (map)",
           "slope S2 (lact)", "slope S3 (lact)", "slope S4 (lact)", "slope S5 (lact)",
           paste0("Upsilon (", 1:16, ", ", rep(1:16, each = 16), ")"), 
           "AR (hemo)", "AR (hr)", "AR (MAP)", "AR (lact)",
           "Var(hemo)", "Cov(hemo, hr)", "Cov(hemo, map)", "Cov(hemo, lact)", 
           "Cov(hr, hemo)", "Var(hr)", "Cov(hr, map)", "Cov(hr, lact)",
           "Cov(map, hemo)", "Cov(map, hr)", "Var(map)", "Cov(map, lact)",
           "Cov(lact, hemo)", "Cov(lact, hr)", "Cov(lact, map)", "Var(lact)",
           "intercept: S1 --> S2", "slope: S1 --> S2", "intercept: S1 --> S4", "slope: S1 --> S4",
           "intercept: S2 --> S3", "slope: S2 --> S3", "intercept: S2 --> S4", "slope: S2 --> S4",
           "intercept: S3 --> S1", "slope: S3 --> S1", "intercept: S3 --> S2", "slope: S3 --> S2", 
           "intercept: S3 --> S4", "slope: S3 --> S4", "intercept: S4 --> S2", "slope: S4 --> S2", 
           "intercept: S4 --> S5", "slope: S4 --> S5", "intercept: S5 --> S1", "slope: S5 --> S1", 
           "intercept: S5 --> S2", "slope: S5 --> S2", "intercept: S5 --> S4", "slope: S5 --> S4", 
           "logit Pr(init S2)", "logit Pr(init S3)",
           "logit Pr(init S4)", "logit Pr(init S5)",
           paste0("mean (hr): ", Dn_omega_names[1:34]), 
           paste0("mean (map): ", Dn_omega_names[35:84]), 
           paste0("log Upsilon (hr): ", Dn_omega_names[1:34]), 
           paste0("log Upsilon (map): ", Dn_omega_names[35:84]),
           "G(hemo)", "G(hemo, hr)", "G(hemo, map)", "G(hemo, lact)", 
           "G(hr, hemo)", "G(hr)", "G(hr, map)", "G(hr, lact)",
           "G(map, hemo)", "G(map, hr)", "G(map)", "G(map, lact)",
           "G(lact, hemo)", "G(lact, hr)", "G(lact, map)", "G(lact)")

if(simulation) {
    # Estimate the initial state probabilities
    init_prob_mat = matrix(nrow = length(index_seeds), ncol = 5)
    c_s = 1
    for(seed in index_seeds) {
        load(paste0('Data/sim_data_', seed, '.rda'))
        first_ind = c(0, which(diff(data_format[,"EID"]) != 0)) + 1
        init_state = data_format[first_ind, "b_true"]
        init_prob_mat[c_s, ] = c(sum(init_state == 1), sum(init_state == 2),
                                 sum(init_state == 3), sum(init_state == 4), 
                                 sum(init_state == 5)) / length(init_state)
        c_s = c_s + 1
    }
    
    init_par_est = colMeans(init_prob_mat)
    true_par[par_index$init[1]] = log(init_par_est[2] / (1 - sum(init_par_est[2:5])))
    true_par[par_index$init[2]] = log(init_par_est[3] / (1 - sum(init_par_est[2:5])))
    true_par[par_index$init[3]] = log(init_par_est[4] / (1 - sum(init_par_est[2:5])))
    true_par[par_index$init[4]] = log(init_par_est[5] / (1 - sum(init_par_est[2:5])))
}

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------

chain_list = list()
par_means = vector(mode = 'list', length = length(index_seeds))
ind = 0
post_means_mat = NULL
covg = NULL

for(seed in index_seeds){
    

    it_seq = 1:(it_num - start_ind + 1)
    covg_val = FALSE
    
    if(simulation) {
        check_name = paste0('Model_out/mcmc_out_',trialNum, '_', seed, 
                            'it', it_num,'_sim.rda')    
    } else {
        check_name = paste0('Model_out/mcmc_out_',trialNum, '_', seed, 
                            'it', it_num,'.rda')
    }
    
    if(file.exists(check_name)) {
        
        for(it in it_seq) {
            
            if(simulation) {
                load(paste0('Model_out/mcmc_out_',trialNum, '_', seed, 
                            'it', it + (start_ind - 1),'_sim.rda'))    
            } else {
                load(paste0('Model_out/mcmc_out_',trialNum, '_', seed, 
                            'it', it + (start_ind - 1),'.rda'))
            }
            
            print(paste0(seed, ": ", it))
            print("accept")
            print(mcmc_out$accept)
            
            par_index = mcmc_out$par_index
            
            if(it == 1) {
                ind = ind + 1
                chain_list[[ind]] = mcmc_out$chain[500:1000, ]
            } else {
                chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out$chain)
            }
            
            rm(mcmc_out)
        }
        
        covg_val = TRUE
        
    } else {
        print(paste0("Missing! ", check_name))
    }
    
    if(covg_val) {
        post_means_temp = matrix(colMeans(chain_list[[ind]]), nrow = 1)
        post_means_mat = rbind(post_means_mat, post_means_temp)
        par_means[[seed]] = c(post_means_temp)
        
        covg_low = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.025)})
        covg_high = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.975)})
        covg_temp = as.numeric(true_par <= covg_high & true_par >= covg_low)
        covg = rbind(covg, covg_temp)    
    }
}

if(!simulation) {save(par_means, file = paste0('Model_out/par_means_it', it_num, '_', as.numeric(simulation), '.rda'))}

upsilon_ind = matrix(1:length(par_index$upsilon), ncol = 16)
upsilon_ind[upper.tri(upsilon_ind, diag = F)] = 0
upsilon_ind = c(upsilon_ind)
upsilon_ind = upsilon_ind[upsilon_ind != 0]


stacked_chains = do.call( rbind, chain_list)

pdf_title = paste0('Plots/trace_plot_', as.numeric(simulation), '_it', it_num, '.pdf')
pdf(pdf_title)
par(mfrow=c(3, 2))
lab_ind = 0
for(s in names(par_index)){
    temp_par = par_index[[s]]
    if(s == "upsilon") {
        temp_par = temp_par[upsilon_ind]
    }
    
    for(r in temp_par){
        # lab_ind = lab_ind + 1
        lab_ind = r
        parMean = round( mean(stacked_chains[,r], na.rm = T), 4)
        parMedian = round( median(stacked_chains[,r], na.rm = T), 4)
        upper = quantile( stacked_chains[,r], prob=.975, na.rm = T)
        lower = quantile( stacked_chains[,r], prob=.025, na.rm = T)
        
        y_limit = range(stacked_chains[,r])
        plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,nrow(chain_list[[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"))
        
        for(seed in 1:length(chain_list)) { lines( chain_list[[seed]][,r], type='l', col=seed) }
        
        if(simulation) {
            x_label = paste0('Mean =',toString(parMean),
                             ' Median =',toString(parMedian),
                             ' True =', round(true_par[r], 3)) 
        } else {
            x_label = paste0('Mean =',toString(parMean),
                             ' Median =',toString(parMedian))     
        }
        
        hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, freq=FALSE,
              xlab=x_label, main = "")
        abline( v=upper, col='red', lwd=2, lty=2)
        abline( v=lower, col='purple', lwd=2, lty=2)
        if(simulation) { abline( v=true_par[r], col='green', lwd=2, lty=2) }
    }   
}

if(simulation) {
    par(mfrow=c(3, 3))
    lab_ind = 0
    for(s in names(par_index)){
        temp_par = par_index[[s]]
        if(s == "upsilon") {
            temp_par = temp_par[upsilon_ind]
        }
        for(r in temp_par){
            lab_ind = r
            boxplot(post_means_mat[,r], main = labels[r],
                    xlab = paste0('95% Covg = ', round(mean(covg[,r]), 4)),
                    ylab = paste0('truth = ', true_par[r]))
            abline(h = true_par[r], col = 'red')
        }
    }
}

dev.off()

# main = paste0("(GR, sGR) = (", round(GR_univ[r,1], digits = 4),
#                                           ", ", round(GR_univ[r,2], digits = 4),"), (",
#                                           round(GR_mult[1], digits = 3),", ",
#                                           round(GR_mult[2], digits = 3),")")



