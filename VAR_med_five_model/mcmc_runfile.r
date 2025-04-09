source('mcmc_routine.r')

# Input will be a number 1-25 (5 seeds with 5 sampling routines)
args = commandArgs(TRUE)

# df_num = as.numeric(args[1])
# seed_num = df_num
# sampling_num = 5
# states_per_step = 0
# steps_per_it = 1

seed_num = as.numeric(args[1])
sampling_num = floor((seed_num - 1) / 5) + 1
seed_num = seed_num - 5 * floor((seed_num - 1)/5)
p = 2
states_per_step = p + 1
steps_per_it = 1
df_num = 1

if(sampling_num %in% c(4,5)) {
    states_per_step = 0
    steps_per_it = p - 1
}

set.seed(seed_num)
steps  = 20000
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
    max_ind = 5
    if(max_ind > 5) burnin = 0
    
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

# Parameter initialization -----------------------------------------------------
load('Data_sim/true_par_index.rda')

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
par[par_index$vec_zeta] = c(-4.7405, 4.5, -5.2152,   1, -3.6473,-0.5, -3.1475, -0.2, 
                            -6.4459,  -1, -3.9404,   2, -4.2151,   1, -4.1778,  2.5, 
                            -3.0523,   0, -6.4459,-0.2, -4.2404, 3.5, -4.2151,    1)

par[par_index$vec_init] = c(-1, 0, -0.5, 0.1)

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
}
# -----------------------------------------------------------------------------
A = list()
W = list()

for(i in 1:length(EIDs)){
    if(simulation) {
        A[[i]] = alpha_i_mat[[i]]
        W[[i]] = omega_i_mat[[i]]
    } else {
        A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol =1)
        W[[i]] = matrix(par[par_index$omega_tilde], ncol =1)
    }
}

B = list()
if(max_ind > 5) {
    
} else {
    for(ii in 1:length(EIDs)) {
        i = EIDs[ii]
        
        # if(simulation) {
        #     # Initialize at the "Maximum likelihood state sequence"
        #     
        #     # Initialize at the "true" state sequence
        #     B[[ii]] = matrix(b_chain[data_format[,"EID"] == i], ncol = 1)
        # } else {
            B[[ii]] = matrix(rep(1, sum(data_format[,"EID"] == i)), ncol = 1)
        # }
    }
}
# -----------------------------------------------------------------------------

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, seed_num, 
                         trialNum, Dn_omega, simulation, bleed_indicator, 
                         max_ind, df_num, sampling_num, states_per_step, steps_per_it)
e_time = Sys.time() - s_time; print(e_time) 