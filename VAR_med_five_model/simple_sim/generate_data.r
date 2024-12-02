set.seed(2024)
N = 100
n_state = 5

#          1->2,    1->4,    2->3,    2->4,    3->1,    3->2,    3->4,    4->2,
#          4->5,    5->1,    5->2,    5->4
par = c(0, 1, -1, 2, -2,
        -4.7405, -5.2152, -3.6473, -3.1475, -6.4459, -3.9404, -4.2151, -4.1778,
        -3.0523, -6.4459, -4.2404, -4.2151)
par_index = list()
par_index$mu = 1:5
par_index$t_p = 6:17

zeta = par[par_index$t_p]
zeta = exp(zeta)

Q = matrix(c(      1,  zeta[1],       0,  zeta[2],       0,
                   0,        1, zeta[3],  zeta[4],       0,
             zeta[5],  zeta[6],       1,  zeta[7],       0,
                   0,  zeta[8],       0,        1, zeta[9],
            zeta[10], zeta[11],       0, zeta[12],       1), ncol=5, byrow=T)

P = Q / rowSums(Q)

init_prob = c(0.245, 0.090, 0.245, 0.149, 0.271)

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
save(data_format, file = 'Data/data_format.rda')
