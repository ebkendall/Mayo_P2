set.seed(2024)
N = 100

# par = c(1, 0, 0.4054651, -0.4054651)
par = c(1, 0, 2.197225, -2.197225)
par_index = list()
par_index$mu = 1:2
par_index$t_p = 3:4

P = matrix(c(                         1, exp(par[par_index$t_p][1]),
             exp(par[par_index$t_p][2]),                         1),
             ncol = 2, byrow = T)

P = P / rowSums(P)

data_format = NULL
for(i in 1:N) {
    n_i = rpois(n = 1, lambda = 1000)
    
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

# #          1->2,    1->4,    2->3,    2->4,    3->1,    3->2,    3->4,    4->2,
# #          4->5,    5->1,    5->2,    5->4
# par = c(0.4054651, 0.4054651,
#         -4.7405, -5.2152, -3.6473, -3.1475, -6.4459, -3.9404, -4.2151, -4.1778, 
#         -3.0523, -6.4459, -4.2404, -4.2151)
# par_index = list()
# par_index$mu = 1:2
# par_index$t_p = 3:14
# 
# zeta = par[par_index$t_p]
# 
# Q = matrix(c(   1,   q1,  0,  q2,  0,
#                 0,    1, q3,  q4,  0,
#                 q5,   q6,  1,  q7,  0,
#                 0,   q8,  0,   1, q9,
#                 q10,  q11,  0, q12,  1), ncol=5, byrow=T)
# 
# P_i = Q / rowSums(Q)
# P = matrix(c(                         1, exp(par[par_index$t_p][1]),
#                                       exp(par[par_index$t_p][1]),                         1),
#            ncol = 2, byrow = T)
# 
# P = P / rowSums(P)
