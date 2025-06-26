index_seeds = c(1:100)
it_num = 1

# Parameter initialization -----------------------------------------------------
par_index = list()
par_index$alpha = 1:8
par_index$zeta = 9:12
par_index$diag_R = 13:16
par_index$init = 17:18
par_index$diag_G = 19:22

true_par = rep(0, tail(par_index$init, 1))
true_par[par_index$alpha] = c( -5,   5,
                               10, -10,
                               -10,  10,
                               5,  -5)
true_par[par_index$zeta] = c(-2, -2, -1.5, -1.5)
true_par[par_index$diag_R] = c(1.386294, 1.386294, 1.386294, 1.386294)
true_par[par_index$init] = c(0, 0)
true_par[par_index$diag_G] = c(1.386294, 1.386294, 1.386294, 1.386294)

labels = c("S2 slope y1", "S3 slope y1",  
           "S2 slope y2", "S3 slope y2",
           "S2 slope y3", "S3 slope y3",
           "S2 slope y4", "S3 slope y4",
           "logit baseline 1 -> 2", "logit baseline 2 -> 3",
           "logit baseline 3 -> 1", "logit baseline 3 -> 2",
           "log R(1,1)", "log R(2,2)", "log R(3,3)", "log R(4,4)",
           "logit init S2", "logit init S3", 
           "log G(1,1)", "log G(2,2)", "log G(3,3)", "log G(4,4)")

# Estimate the initial state probabilities
init_prob_mat = matrix(nrow = length(index_seeds), ncol = 3)
for(seed in index_seeds) {
    load(paste0('Data/data_format', seed, '.rda'))
    first_ind = c(0, which(diff(data_format[,"id"]) != 0)) + 2 # ignore the first state
    init_state = data_format[first_ind, "state"]
    init_prob_mat[seed, ] = c(sum(init_state == 1), sum(init_state == 2), sum(init_state == 3)) / length(init_state)
}

init_par_est = colMeans(init_prob_mat)
true_par[par_index$init[1]] = log(init_par_est[2] / (1 - init_par_est[2] - init_par_est[3]))
true_par[par_index$init[2]] = log(init_par_est[3] / (1 - init_par_est[2] - init_par_est[3]))

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))

ind = 0
post_means_mat = NULL
covg = NULL

for(seed in index_seeds){
    
    it_seq = 1:it_num
    covg_val = FALSE
    
    for(it in it_seq) {
        file_name = paste0('Model_out/mcmc_out_',seed, '_it_', it,'.rda')
        if(file.exists(file_name)) {
            load(file_name)
            print(paste0(seed, ": ", file_name))
            print("accept")
            print(mcmc_out$accept)
            
            par_index = mcmc_out$par_index
            
            if(it == 1) {
                ind = ind + 1
                chain_list[[ind]] = mcmc_out$chain[1:1000, ]
            } else {
                chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out$chain)
            }
            
            rm(mcmc_out)
            
            covg_val = TRUE
        } else {
            print(paste0("Missing! ", file_name))
        }
        
    }
    
    if(covg_val) {
        post_means_temp = matrix(colMeans(chain_list[[ind]]), nrow = 1)
        post_means_mat = rbind(post_means_mat, post_means_temp)
        
        covg_low = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.025)})
        covg_high = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.975)})
        covg_temp = as.numeric(true_par <= covg_high & true_par >= covg_low)
        covg = rbind(covg, covg_temp)    
    }
}

# # Compute the OG Gelman-Rubin stat. for testing convergence --------------------
# # Univariate --------------------
# gr_stat = rep(0, length(true_par))
# n = nrow(chain_list[[1]])
# m = length(index_seeds)
# for(p in 1:length(true_par)) {
#     x_i_bar = s_i = rep(0, length(index_seeds))
#     for(i in 1:m) {
#         x_i_bar[i] = mean(chain_list[[i]][,p])
#         s_i[i] = (1/(n - 1)) * sum((chain_list[[i]][,p] - x_i_bar[i])^2)
#     }
#     mu_hat = mean(x_i_bar)
#     s_2 = mean(s_i)
    
#     B_n = (1 / (m - 1)) * sum((x_i_bar - mu_hat)^2)
#     sigma_hat_2 = ((n-1)/n) * s_2 + B_n
#     gr_stat[p] = sqrt(sigma_hat_2 / s_2)
# }

# # Multivariate ------------------
# out_prod = function(x) {
#     c1 = matrix(x, ncol = 1)
#     return(c1 %*% t(c1))
# }
# X_bar_means = matrix(nrow = length(true_par), ncol = m)
# S_i_list = vector(mode = 'list', length = m)
# for(i in 1:m) {
#     X_bar_means[,i] = colMeans(chain_list[[i]])
#     X_diff = sweep(chain_list[[i]], 2, colMeans(chain_list[[i]]), '-')
#     X_diff_df = as.data.frame(t(X_diff))
    
#     X_outer = lapply(X_diff_df, out_prod)
#     S_i_list[[i]] = (1 / (n-1)) * Reduce('+', X_outer)
# }
# big_S = (1 / m) * Reduce('+', S_i_list)
# vec_mu_hat = apply(X_bar_means, 1, mean)

# B_diff = sweep(X_bar_means, 1, vec_mu_hat, '-')
# B_diff_df = as.data.frame(B_diff)
# B_outer = lapply(B_diff_df, out_prod)
# big_B_n = (1 / (m-1)) * Reduce('+', B_outer)

# lambda_max = abs(eigen(solve(big_S) %*% (n * big_B_n))$values[1])
# R_hat = sqrt(((n-1) / n) + (lambda_max / n))

# # Compute the stable Gelman-Rubin stat. for testing convergence ----------------
# library(stableGR)
# R_hat_stable = stable.GR(x = chain_list,
#                          multivariate = T,
#                          mapping = "determinant",
#                          method = "lug",
#                          size = NULL,
#                          autoburnin = FALSE,
#                          blather = FALSE)

# GR_univ = cbind(gr_stat, R_hat_stable$psrf)
# GR_mult = c(R_hat, R_hat_stable$mpsrf)


stacked_chains = do.call( rbind, chain_list)

pdf_title = paste0('Plots/trace_plot.pdf')
pdf(pdf_title)
par(mfrow=c(3, 2))
lab_ind = 0
for(s in names(par_index)){
    temp_par = par_index[[s]]
    for(r in temp_par){
        # lab_ind = lab_ind + 1
        lab_ind = r
        parMean = round( mean(stacked_chains[,r], na.rm = T), 4)
        parMedian = round( median(stacked_chains[,r], na.rm = T), 4)
        upper = quantile( stacked_chains[,r], prob=.975, na.rm = T)
        lower = quantile( stacked_chains[,r], prob=.025, na.rm = T)
        
        title_color = "black"
        
        y_limit = range(stacked_chains[,r])
        plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,nrow(chain_list[[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"),
              col.main = title_color)
        
        for(seed in 1:length(chain_list)) { lines( chain_list[[seed]][,r], type='l', col=seed) }
        
        x_label = paste0('Mean =',toString(parMean),
                         ' Median =',toString(parMedian),
                         ' True =', round(true_par[r], 3)) 
        
        hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, freq=FALSE,
              xlab=x_label, main = "")
        abline( v=upper, col='red', lwd=2, lty=2)
        abline( v=lower, col='purple', lwd=2, lty=2)
        abline( v=true_par[r], col='green', lwd=2, lty=2)
    }   
}

par(mfrow=c(3, 3))
lab_ind = 0
for(s in names(par_index)){
    temp_par = par_index[[s]]
    for(r in temp_par){
        lab_ind = r
        boxplot(post_means_mat[,r], main = labels[r],
                xlab = paste0('95% Covg = ', round(mean(covg[,r]), 4)),
                ylab = paste0('truth = ', true_par[r]))
        abline(h = true_par[r], col = 'red')
    }
}

dev.off()

# main = paste0("(GR, sGR) = (", round(GR_univ[r,1], digits = 4),
#                                           ", ", round(GR_univ[r,2], digits = 4),"), (",
#                                           round(GR_mult[1], digits = 3),", ",
#                                           round(GR_mult[2], digits = 3),")")



