index_seeds = c(1:10)
trialNum = 1
simulation = F

it_num = 2
start_ind = 1

# Parameter initialization -----------------------------------------------------
par_index = list()
par_index$beta = 1:4
par_index$alpha_tilde = 5:20
par_index$upsilon = 21:276
par_index$A = 277:296
par_index$R = 297:312
par_index$zeta = 313:336
par_index$init = 337:340
par_index$omega_tilde = 341:424
par_index$G = 425:440

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
true_par[par_index$A] = c(rep(2, 4), rep(-2, 4), rep(0, 4), rep(-2, 4), rep(0, 4))
true_par[par_index$R] = c(diag(c(9, 9, 9, 9)))
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
true_par[par_index$G] = c(diag(c(9, 9, 9, 9)))

load('Data/Dn_omega_names.rda')

labels = c("beta (hemo)", "beta (hr)", "beta (map)", "beta (lact)",
           "slope S2 (hemo)", "slope S3 (hemo)", "slope S4 (hemo)", "slope S5 (hemo)",
           "slope S2 (hr)", "slope S3 (hr)", "slope S4 (hr)", "slope S5 (hr)",
           "slope S2 (map)", "slope S3 (map)", "slope S4 (map)", "slope S5 (map)",
           "slope S2 (lact)", "slope S3 (lact)", "slope S4 (lact)", "slope S5 (lact)",
           paste0("Upsilon (", 1:16, ", ", rep(1:16, each = 16), ")"), 
           "AR (hemo, S1)", "AR (hr, S1)", "AR (map, S1)", "AR (lact, S1)",
           "AR (hemo, S2)", "AR (hr, S2)", "AR (map, S2)", "AR (lact, S2)",
           "AR (hemo, S3)", "AR (hr, S3)", "AR (map, S3)", "AR (lact, S3)",
           "AR (hemo, S4)", "AR (hr, S4)", "AR (map, S4)", "AR (lact, S4)",
           "AR (hemo, S5)", "AR (hr, S5)", "AR (map, S5)", "AR (lact, S5)",
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
A_list = list()
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
                A_list[[ind]] = t(do.call(cbind, mcmc_out$alpha_i))
            } else {
                chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out$chain)
                A_list[[ind]] = rbind(A_list[[ind]], t(do.call(cbind, mcmc_out$alpha_i)))
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
        
        print(paste0(labels[lab_ind], ": ", round(parMedian, 3), " [", round(lower, 3), ", ", round(upper, 3), "]"))
        
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

# Plot the sampled alpha_i
alpha_i_names = c("slope S2 (hemo)", "slope S3 (hemo)", "slope S4 (hemo)", "slope S5 (hemo)",
                  "slope S2 (hr)", "slope S3 (hr)", "slope S4 (hr)", "slope S5 (hr)",
                  "slope S2 (map)", "slope S3 (map)", "slope S4 (map)", "slope S5 (map)",
                  "slope S2 (lact)", "slope S3 (lact)", "slope S4 (lact)", "slope S5 (lact)")
par(mfrow = c(2,2))
for(a in 1:16) {
    all_a = NULL
    freq_counts = NULL
    for(i in 1:length(A_list)) {
        all_a = cbind(all_a, A_list[[i]][,a])
        
        temp_hist = hist(A_list[[i]][,a], breaks = sqrt(nrow(A_list[[i]])), plot = F)
        freq_counts = c(freq_counts, max(temp_hist$counts))
    }
    
    x_margin = c(min(c(all_a)), max(c(all_a)))
    y_margin = c(0, max(freq_counts))
    
    hist(all_a[,1], main = alpha_i_names[a], xlab = paste0('alpha_i(', a, ')'),
         breaks = sqrt(nrow(all_a)), col = 1, xlim = x_margin, ylim = y_margin)
    for(i in 2:length(A_list)) {
        hist(all_a[,i], breaks = sqrt(nrow(all_a)), col = i, add = T)
    }
}

dev.off()

# main = paste0("(GR, sGR) = (", round(GR_univ[r,1], digits = 4),
#                                           ", ", round(GR_univ[r,2], digits = 4),"), (",
#                                           round(GR_mult[1], digits = 3),", ",
#                                           round(GR_mult[2], digits = 3),")")



