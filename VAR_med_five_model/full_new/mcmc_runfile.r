source('mcmc_routine.r')

args = commandArgs(TRUE)
seed_num = as.numeric(args[1])
set.seed(seed_num)

simulation = T

# Load data --------------------------------------------------------------------
data_format = NULL

if(simulation) {
    trialNum = 1
    max_ind = 5
    
    load(paste0('Data/sim_data_', seed_num, '.rda'))
    print(paste0('SIM: seed ', seed_num, ' trial ', trialNum))
} else {
    trialNum = 1
    max_ind = 5
    
    load('Data/data_format_train_update.rda')
    print(paste0('REAL: seed ', seed_num, ' trial ', trialNum))
}

Y = data_format[, c('EID', 'hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
x = data_format[,c('n_RBC_admin'), drop=F]
z = cbind(1, data_format[,c('RBC_ordered'), drop=F])

EIDs = unique(data_format[,'EID'])

# Parameter initialization -----------------------------------------------------
par_index = list()
par_index$beta = 1:4
par_index$alpha_tilde = 5:20
par_index$upsilon = 21:276
par_index$A = 277:296
par_index$R = 297:312
par_index$zeta = 313:336
par_index$init = 337:340
par_index$omega_tilde = 341:424
par_index$G = 425:440

par = rep(0, max(do.call('c', par_index)))

if(simulation) {
    load('Model_out/mcmc_out_1_1it2.rda')
    par[par_index$beta] = mcmc_out$chain[nrow(mcmc_out$chain), mcmc_out$par_index$beta]
    par[par_index$alpha_tilde] = mcmc_out$chain[nrow(mcmc_out$chain), mcmc_out$par_index$alpha_tilde]
    par[par_index$upsilon] = mcmc_out$chain[nrow(mcmc_out$chain), mcmc_out$par_index$upsilon]
    par[par_index$A] = mcmc_out$chain[nrow(mcmc_out$chain), mcmc_out$par_index$A]
    par[par_index$R] = mcmc_out$chain[nrow(mcmc_out$chain), mcmc_out$par_index$R]
    par[par_index$zeta] = mcmc_out$chain[nrow(mcmc_out$chain), mcmc_out$par_index$zeta]
    par[par_index$init] = mcmc_out$chain[nrow(mcmc_out$chain), mcmc_out$par_index$init]
    par[par_index$omega_tilde] = mcmc_out$chain[nrow(mcmc_out$chain), mcmc_out$par_index$omega_tilde]
    par[par_index$G] = c(diag(c(9, 9, 9, 9)))
    
    init_par_est = c(312,  25,  27,  69,  67); init_par_est = init_par_est/sum(init_par_est)
    par[par_index$init][1] = log(init_par_est[2] / (1 - sum(init_par_est[2:5])))
    par[par_index$init][2] = log(init_par_est[3] / (1 - sum(init_par_est[2:5])))
    par[par_index$init][3] = log(init_par_est[4] / (1 - sum(init_par_est[2:5])))
    par[par_index$init][4] = log(init_par_est[5] / (1 - sum(init_par_est[2:5])))
    
    rm(mcmc_out)
} else {
    par[par_index$beta] = c(0.25, -2, 2, -0.25) 
    par[par_index$alpha_tilde] = c(-1,  1, 0, 0,
                                    7, -7, 0, 0,
                                   -7,  7, 0, 0,
                                    1, -1, 0, 0)
    par[par_index$upsilon] = c(diag(c(0.25, 0.25,  4,  4,
                                      2.25, 2.25, 25, 25,
                                      2.25, 2.25, 25, 25,
                                      0.25, 0.25,  4,  4)))
    par[par_index$A] = c(rep(2, 4), rep(-2, 4), rep(0, 4), rep(-2, 4), rep(0, 4))
    par[par_index$R] = c(diag(c(0.5, 1.5, 1.5, 0.5)))
    #    transitions:          1->2,         1->4,         2->3,         2->4, 
    #                          3->1,         3->2,         3->4,         4->2, 
    #                          4->5,         5->1,         5->2,         5->4
    par[par_index$zeta] = c(-7.2405, 2.5, -6.2152,   1, -2.6473,  -1, -6.1475,  -1, 
                            -9.4459,  -1, -7.2404, 2.5, -7.2151,   1, -7.1778, 2.5, 
                            -5.2151,   0, -9.4459,  -1, -7.2404, 2.5, -5.2151,   0)
    par[par_index$init] = c(0, 0, 0, 0)
    par[par_index$omega_tilde]= 2 * c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                                      -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                                      -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                                      -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
    par[par_index$G] = c(diag(c(8, 32, 32, 8)))
}

n_state = 5

A = list()
B = list()

if(simulation) {
    load(paste0('Data/alpha_i_mat_', seed_num, '.rda'))
    load('Data/Dn_omega_sim.rda')
    load(paste0('Data/bleed_indicator_sim_', seed_num, '.rda'))
    
    A = alpha_i_mat
    Dn_omega = Dn_omega_sim
    b_chain = data_format[, "b_true"]
    
    for(ii in 1:length(EIDs)){
        i = EIDs[ii]
        
        B[[ii]] = matrix(b_chain[data_format[,"EID"] == i], ncol = 1)
    }
} else {
    load('Data/Dn_omega_update.rda')
    bleed_indicator = b_ind_fnc(data_format)

    if(max_ind > 5) {

        # 9: from it 1-3
        # all seeds for rest
        chosen_seed = 3
        load(paste0('Model_out/mcmc_out_', trialNum, '_', chosen_seed, 'it',
                    max_ind - 5, '.rda'))

        par = mcmc_out$chain[nrow(mcmc_out$chain), ]
        par[par_index$A] = c(rep(2, 4), rep(-2, 4), rep(0, 4), rep(-2, 4), rep(0, 4))
        par[par_index$R] = c(diag(c(0.5, 1.5, 1.5, 0.5)))
        par[par_index$upsilon] = c(diag(c(0.25, 0.25,  4,  4,
                                          2.25, 2.25, 25, 25,
                                          2.25, 2.25, 25, 25,
                                          0.25, 0.25,  4,  4)))
        par[par_index$omega_tilde]= 2 * c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                                          -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                                          -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                                          -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
        
        b_chain = mcmc_out$B_chain[nrow(mcmc_out$B_chain), ]

        A = mcmc_out$alpha_i

        for(ii in 1:length(EIDs)){
            i = EIDs[ii]

            B[[ii]] = matrix(b_chain[data_format[,"EID"] == i], ncol = 1)
        }

        rm(mcmc_out)
    } else {
        
        b_chain = rep(1, nrow(data_format))

        for(ii in 1:length(EIDs)){
            i = EIDs[ii]

            A[[ii]] = matrix(par[par_index$alpha_tilde], ncol =1)
            B[[ii]] = matrix(b_chain[data_format[,"EID"] == i], ncol = 1)
        }
    }
}
# -----------------------------------------------------------------------------

steps  = 50000
burnin =  5000

# if(max_ind > 5) {burnin = 0}

s_time = Sys.time()
mcmc_out = mcmc_routine(steps, burnin, seed_num, trialNum, simulation, max_ind,
                        par, par_index, Y, x, z, B, A, Dn_omega, bleed_indicator)
e_time = Sys.time() - s_time; print(e_time) 
