source('mcmc_routine.r')

# Input will be a number 1-25 (5 seeds with 5 sampling routines)
args = commandArgs(TRUE)

df_num = as.numeric(args[1])
seed_num = df_num
sampling_num = 5
before_t1 = TRUE # Do we handle the state changes before t1 or not?

# seed_num = as.numeric(args[1])
# sampling_num = floor((seed_num - 1) / 5) + 1
# seed_num = seed_num - 5 * floor((seed_num - 1)/5)
# p = 2
# states_per_step = p + 1
# steps_per_it = 1
# df_num = 1
# 
# if(sampling_num %in% c(4,5)) {
#     states_per_step = 0
#     steps_per_it = p - 1
# }

set.seed(seed_num)
steps  = 10000
burnin =  5000

simulation = T
data_format = NULL

if(simulation) {
    trialNum = 1
    max_ind = 5
    
    load(paste0('Data_sim/use_data_', df_num, '.rda'))
    data_format = use_data
    
    print(paste0('SIM: seed ', seed_num, ' samp ', 
                 sampling_num, ' trial ', trialNum))
} else {
    trialNum = 1
    max_ind = 10
    if(max_ind > 5) {burnin = 0}
    
    load('Data_real/data_format_train.rda')
    print(paste0('REAL: seed ', seed_num, ' samp ', 
                 sampling_num, ' trial ', trialNum))
}

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = unique(data_format[,'EID'])

# Covariate on the mean process
x = data_format[,c('n_RBC_admin'), drop=F]

# Covariate on the state process
z = cbind(1, data_format[,c('RBC_ordered'), drop=F])

# Indexing initialization ------------------------------------------------------
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:24
par_index$vec_sigma_upsilon = 25:424
par_index$vec_A = 425:428
par_index$vec_R = 429:444
par_index$vec_zeta = 445:468
par_index$vec_init = 469:472
par_index$omega_tilde = 473:556
par_index$vec_upsilon_omega = 557:640

# Parameter initialization -----------------------------------------------------
par = rep(0, max(par_index$vec_upsilon_omega))

par[par_index$vec_beta] = c(0.25, -2, 2, -0.25) # one unit of RBC -> 1 unit increase in hemo in 1 hour
par[par_index$vec_alpha_tilde] = c( 9, -1,  1, 0, 0,
                                   85,  5, -5, 0, 0,
                                   75, -5,  5, 0, 0,
                                    5,  1, -1, 0, 0)
par[par_index$vec_sigma_upsilon] = c(diag(c(  4, 0.01, 0.01, 0.25, 0.25, 
                                              100,    1,    1,   25,   25, 
                                              100,    1,    1,   25,   25, 
                                              1, 0.01, 0.01, 0.25, 0.25)))
par[par_index$vec_A] = c(rep(0, 4))
par[par_index$vec_R] = c(diag(c(9, 9, 9, 9)))

#    transitions:              1->2,         1->4,         2->3,         2->4, 
#                              3->1,         3->2,         3->4,         4->2, 
#                              4->5,         5->1,         5->2,         5->4,
par[par_index$vec_zeta] = c(-3.7405, 2.5, -4.2152,   1, -2.6473,-0.5, -2.1475, -0.2, 
                            -3.4459,  -1, -2.9404,   1, -3.2151,   1, -3.1778,  1.5, 
                            -2.0523,   0, -3.4459,-0.2, -3.2404, 2.5, -3.2151,    1)

par[par_index$vec_init] = c(-5, -5, -5, -5)

par[par_index$omega_tilde]= 2 * c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                                  -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                                  -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                                  -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)

par[par_index$vec_upsilon_omega] = rep(1, length(par_index$vec_upsilon_omega))
# -----------------------------------------------------------------------------
if(simulation) {
    load('Data_sim/true_pars.rda')
    load(paste0('Data_sim/alpha_i_mat_', df_num, '.rda'))
    load(paste0('Data_sim/omega_i_mat_', df_num, '.rda'))
    load('Data_sim/Dn_omega_sim.rda')
    load(paste0('Data_sim/bleed_indicator_sim_', df_num, '.rda'))
    
    par = true_pars
    Dn_omega = Dn_omega_sim
    
    b_chain = data_format[, "b_true"]
} else {
    load('Data_real/Dn_omega.rda')
    bleed_indicator = b_ind_fnc(data_format)
    
    if(max_ind > 5) {
        # it5, samp 4: seed 4
        # it5, samp 5: seed 1
        if(sampling_num == 4) {
            chosen_seed = 4
        } else if(sampling_num == 5) {
            chosen_seed = 1
        }
        
        load(paste0('Model_out/mcmc_out_', trialNum, '_', chosen_seed, 'it', 
                    max_ind - 5, '_samp', sampling_num, '.rda'))
        par = mcmc_out_temp$chain[nrow(mcmc_out_temp$chain), ]
        b_chain = mcmc_out_temp$B_chain[nrow(mcmc_out_temp$B_chain), ]
        
        rm(mcmc_out_temp)
    } else {
        b_chain = rep(1, nrow(data_format))
    }
}
# -----------------------------------------------------------------------------
A = list()
W = list()
B = list()

for(ii in 1:length(EIDs)){
    i = EIDs[ii]
    
    if(simulation) {
        A[[ii]] = alpha_i_mat[[ii]]
        W[[ii]] = omega_i_mat[[ii]]
    } else {
        A[[ii]] = matrix(par[par_index$vec_alpha_tilde], ncol =1)
        W[[ii]] = matrix(par[par_index$omega_tilde], ncol =1)
    }
    
    B[[ii]] = matrix(b_chain[data_format[,"EID"] == i], ncol = 1)
}
# -----------------------------------------------------------------------------

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, seed_num, 
                         trialNum, Dn_omega, simulation, bleed_indicator, 
                         max_ind, df_num, sampling_num, before_t1)
e_time = Sys.time() - s_time; print(e_time) 