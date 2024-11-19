source('mcmc_routine.r')

args = commandArgs(TRUE)
seed_num = as.numeric(args[1])
sampling_num = NULL

if(seed_num <= 3) {
    sampling_num = 1
} else if(seed_num > 3 & seed_num <= 6) {
    seed_num = seed_num - 3
    sampling_num = 2
} else if(seed_num > 6 & seed_num <= 9) {
    seed_num = seed_num - 6
    sampling_num = 3
} else {
    seed_num = seed_num - 9
    sampling_num = 4
}

set.seed(seed_num)
ind = seed_num

load('Data/data_format.rda')
y = data_format[,"y"]
ids = data_format[,"id"]
EIDs = unique(data_format[,"id"])


# Parameter initialization -----------------------------------------------------
par = c(2, -2, 0.4054651, 0.4054651)
par_index = list()
par_index$mu = 1:2
par_index$t_p = 3:4
# -----------------------------------------------------------------------------

B = list()
P = matrix(c(1, exp(par[par_index$t_p][1]), exp(par[par_index$t_p][1]), 1),
           ncol = 2, byrow = T)

P = P / rowSums(P)

for(i in EIDs){
    
    n_i = sum(data_format[,"id"] == i)
    b_temp = sample(1:2, size = 1, prob = c(0.5, 0.5))
    
    for(k in 2:n_i) {
        b_temp = c(b_temp, sample(1:2, size=1, prob=P[tail(b_temp,1),]))
    }
    
    B[[i]] = matrix(b_temp, ncol = 1)
}
# -----------------------------------------------------------------------------

steps = 20000
burnin = 5000

s_time = Sys.time()
mcmc_out = mcmc_routine(par, par_index, B, y, ids, steps, burnin, ind, sampling_num)
e_time = Sys.time() - s_time; print(e_time)