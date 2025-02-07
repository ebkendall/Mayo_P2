args = commandArgs(TRUE)
seed_num = as.numeric(args[1])

# seed_num = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

sampling_num = NULL
pseudo = F
EM = F

if(pseudo) {
    source('mcmc_routine_pseudo.r')
} else if(EM) {
    source('mcmc_routine_EM.r')
} else {
    source('mcmc_routine.r')
    
    sampling_num = floor((seed_num - 1) / 3) + 1
    seed_num = seed_num - 3 * floor((seed_num - 1)/3)
}

# Num. states sampled per step, Num. steps per MCMC it -------------------------
for(states_per_step in 1:3) {
    # states_per_step = 1
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

    par = c(0, 1, -3, -3)
    par_index = list()
    par_index$mu = 1:2
    par_index$t_p = 3:4
    
    # -----------------------------------------------------------------------------
    
    n_state = 2
    zeta = par[par_index$t_p]
    zeta = exp(zeta)
    Q = matrix(c(      1,  zeta[1],
                 zeta[2],        1), ncol=2, byrow=T)
    
    P = Q / rowSums(Q)
    
    init_prob = rep(1, n_state);
    init_prob = init_prob / sum(init_prob)
    
    B = list()
    for(i in EIDs){
        B[[i]] = matrix(data_format[data_format[,"id"] == i,"state"], ncol = 1)
    }
    # -----------------------------------------------------------------------------
    
    steps  = 30000
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
    
    if(sampling_num == 4) { break }
}