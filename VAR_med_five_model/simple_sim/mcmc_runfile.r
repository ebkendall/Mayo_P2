seed_num = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
sampling_num = NULL
pseudo = F
EM = F

if(pseudo) {
    source('mcmc_routine_pseudo.r')
} else if(EM) {
    source('mcmc_routine_EM.r')
} else {
    source('mcmc_routine.r')
    
    if(seed_num <= 3) {
        sampling_num = 1
    } else if(seed_num > 3 & seed_num <= 6) {
        seed_num = seed_num - 3
        sampling_num = 2
    } else if(seed_num > 6 & seed_num <= 9) {
        seed_num = seed_num - 6
        sampling_num = 3
    } else if(seed_num > 9 & seed_num <= 12) {
        seed_num = seed_num - 9
        sampling_num = 4
    } else {
        seed_num = seed_num - 12
        sampling_num = 5
    }
}

# Num. states sampled per step, Num. steps per MCMC it -------------------------
states_per_step = 2
steps_per_it = 1

set.seed(seed_num)
ind = seed_num

if(EM) {
    load(paste0('Data/data_format', seed_num, '.rda'))
} else {
    load('Data/data_format1.rda')
}
y = data_format[,"y"]
ids = data_format[,"id"]
EIDs = unique(data_format[,"id"])

# Parameter initialization ----------------------------------------------------
par = c(0.5, 0,
        0.405, -0.405)
par_index = list()
par_index$mu = 1:2
par_index$t_p = 3:4

# -----------------------------------------------------------------------------

n_state = 2
zeta = par[par_index$t_p]
zeta = exp(zeta)
Q = matrix(c(       1,  zeta[1],
                    zeta[2],        1), ncol=2, byrow=T)

P = Q / rowSums(Q)

init_prob = c(0.5,0.5)

B = list()
for(i in EIDs){
    
    # n_i = sum(data_format[,"id"] == i)
    # b_temp = sample(1:n_state, size = 1, prob = init_prob)
    # 
    # for(k in 2:n_i) {
    #     b_temp = c(b_temp, sample(1:n_state, size=1, prob=P[tail(b_temp,1),]))
    # }
    # 
    # B[[i]] = matrix(b_temp, ncol = 1)
    B[[i]] = matrix(data_format[data_format[,"id"] == i,"state"], ncol = 1)
}
max_log_posterior = log_post_cpp(as.numeric(EIDs), par, par_index, B, y, ids, 10)
print("Maximum log posterior: "); print(max_log_posterior)
# -----------------------------------------------------------------------------

steps = 30000
burnin =  5000

s_time = Sys.time()
if(pseudo) {
    mcmc_out = mcmc_routine_pseudo(par, par_index, y, ids, steps, burnin, 
                                   ind, sampling_num)
} else if(EM) {
    mcmc_out = mcmc_routine_EM(par, par_index, B, y, ids, steps, burnin, ind)
    save(mcmc_out, file = paste0('Model_out/mcmc_out_',ind,'_EM.rda'))
} else {
    mcmc_out = mcmc_routine(par, par_index, B, y, ids, steps, burnin, ind, 
                            sampling_num, states_per_step, steps_per_it)
}
e_time = Sys.time() - s_time; print(e_time) 