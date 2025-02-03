for(seed_num in 1:100) {
    
    print(seed_num)
    set.seed(seed_num)
    N = 100
    n_state = 2
    
    par = c(0, 1, -3, -3)
    par_index = list()
    par_index$mu = 1:2
    par_index$t_p = 3:4
    
    mu = matrix(par[par_index$mu], ncol = 1)
    zeta = par[par_index$t_p]
    zeta = exp(zeta)
    
    Q = matrix(c(      1,  zeta[1],
                 zeta[2],        1), ncol=2, byrow=T)
    
    P = Q / rowSums(Q)
    
    init_prob = rep(1, n_state);
    init_prob = init_prob / sum(init_prob)
    
    data_format = NULL
    for(i in 1:N) {
        n_i = rpois(n = 1, lambda = 500)
        b_i = NULL
        
        for(k in 1:n_i) {
            if(k == 1) {
                # Initial probabilities are equal
                b_i = sample(1:n_state, size = 1, prob = init_prob)
            } else {
                b_i = c(b_i, sample(1:n_state, size=1, prob=P[tail(b_i,1),]))
            }
        }
        
        D_i = cbind(1, cumsum(b_i == 2))
        mean_i = D_i %*% mu
        y_i = rnorm(n = n_i, mean = mean_i, sd = rep(1,n_i))
        
        sub_i = cbind(i, y_i, b_i)
        data_format = rbind(data_format, sub_i)
    }
    
    colnames(data_format) = c("id", "y", "state")
    save(data_format, file = paste0('Data/data_format', seed_num, '.rda'))    
}
