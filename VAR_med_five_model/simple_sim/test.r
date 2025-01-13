source('mcmc_routine_EM.r')

adjacency_mat = matrix(c(1,1,
                         1,1), ncol = 2, byrow = T)
initialize_cpp(adjacency_mat, 1)

set.seed(2025)
N = 100
n_state = 2

par = c(0.5, 0,
        -1, -4)
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

y = data_format[,"y"]
ids = data_format[,"id"]
EIDs = unique(data_format[,"id"])


B = viterbi_alg(as.numeric(EIDs), par, par_index, y, ids, 10)
B_mle = mle_state_seq(as.numeric(EIDs), par, par_index, y, ids, 10)
b_list = do.call('c', B)
b_list_mle = do.call('c', B_mle)

mean(b_list == data_format[,"state"])
mean(b_list_mle == data_format[,"state"])
