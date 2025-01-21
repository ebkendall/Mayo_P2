N_samp = 100

true_par = c(0.5, 0,
             -1, -4)
par_index = list()
par_index$mu = 1:2
par_index$t_p = 3:4

labels = c("mu1", "mu2", "logit baseline 1 -> 2", "logit baseline 2 -> 1") 

par_chain = matrix(nrow = N_samp, ncol = length(true_par))
state_accuracy = rep(0, N_samp)
state_counts = matrix(ncol=4,nrow = N_samp)
colnames(state_counts) = c('learn1', 'learn2', 'truth1', 'truth2')
for(i in 1:N_samp) {
    load(paste0('Model_out/mcmc_out_', i, '_EM.rda'))
    par_chain[i, ] = tail(mcmc_out$chain, 1)
    b_chain_i = tail(mcmc_out$B_chain, 1)
    
    load(paste0('Data/data_format', i, '.rda'))
    state_accuracy[i] = mean(b_chain_i == data_format[,"state"])
    state_counts[i,] = c(sum(b_chain_i == 1), sum(b_chain_i == 2),
                         sum(data_format[,"state"] == 1), sum(data_format[,"state"] == 2))
}

par(mfrow=c(2,2))
for(i in 1:4) {
    boxplot(par_chain[,i], main = labels[i], xlab = paste0('true = ', true_par[i]))
    abline(h = true_par[i], col = 'green')
}
print(summary(state_accuracy))
print(state_counts)
