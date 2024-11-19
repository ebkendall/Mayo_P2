set.seed(2024)
N = 100

par = c(2, -2, 0.4054651, 0.4054651)
par_index = list()
par_index$mu = 1:2
par_index$t_p = 3:4

P = matrix(c(                         1, exp(par[par_index$t_p][1]),
             exp(par[par_index$t_p][1]),                         1),
             ncol = 2, byrow = T)

P = P / rowSums(P)

data_format = NULL
for(i in 1:N) {
    n_i = rpois(n = 1, lambda = 500)
    
    # Initial probabilities are 50/50
    b_i = sample(1:2, size = 1, prob = c(0.5, 0.5))
    
    for(k in 2:n_i) {
        b_i = c(b_i, sample(1:2, size=1, prob=P[tail(b_i,1),]))
    }
    
    mean_i = par[par_index$mu[b_i]]
    y_i = rnorm(n = n_i, mean = mean_i, sd = rep(1, n_i))
    
    sub_i = cbind(i, y_i, b_i)
    data_format = rbind(data_format, sub_i)
}

colnames(data_format) = c("id", "y", "state")
save(data_format, file = 'Data/data_format.rda')