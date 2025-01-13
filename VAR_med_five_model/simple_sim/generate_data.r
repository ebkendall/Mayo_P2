for(seed_num in 1:100) {
    
    print(seed_num)
    set.seed(seed_num)
    N = 100
    n_state = 2
    
    par = c(1, -1,
            0.405, -0.405)
    par_index = list()
    par_index$mu = 1:2
    par_index$t_p = 3:4
    
    zeta = par[par_index$t_p]
    zeta = exp(zeta)
    
    Q = matrix(c(      1,  zeta[1],
                       zeta[2],        1), ncol=2, byrow=T)
    
    P = Q / rowSums(Q)
    
    init_prob = c(0.5, 0.5)
    
    data_format = NULL
    for(i in 1:N) {
        n_i = rpois(n = 1, lambda = 500)
        
        # Initial probabilities are 50/50
        b_i = sample(1:n_state, size = 1, prob = init_prob)
        
        for(k in 2:n_i) {
            b_i = c(b_i, sample(1:n_state, size=1, prob=P[tail(b_i,1),]))
        }
        
        mean_i = par[par_index$mu[b_i]]
        y_i = rnorm(n = n_i, mean = mean_i, sd = rep(1, n_i))
        
        sub_i = cbind(i, y_i, b_i)
        data_format = rbind(data_format, sub_i)
    }
    
    colnames(data_format) = c("id", "y", "state")
    save(data_format, file = paste0('Data/data_format', seed_num, '.rda'))    
}
