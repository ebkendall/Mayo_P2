source('test_routine.r')

args = commandArgs(TRUE)
seed_num = as.numeric(args[1])
set.seed(seed_num)

# Load data --------------------------------------------------------------------
load('Data/data_format_TEST.rda')

# Remove clinic and RBC rule to leave this as uninformed
data_format[,'RBC_rule'] = 0
data_format[,'clinic_rule'] = 0

load('Data/Dn_omega_TEST.rda')

Y = data_format[, c('EID', 'hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
x = data_format[,c('n_RBC_admin'), drop=F]
z = cbind(1, data_format[,c('RBC_ordered'), drop=F])

EIDs = unique(data_format[,'EID'])

bleed_indicator = b_ind_fnc(data_format)

# Parameter initialization -----------------------------------------------------
par_index = list()
par_index$beta = 1:4
par_index$alpha_tilde = 5:20
par_index$upsilon = 21:276
par_index$A = 277:280
par_index$R = 281:296
par_index$zeta = 297:320
par_index$init = 321:324
par_index$omega_tilde = 325:408
par_index$eta_omega = 409:492
par_index$G = 493:508

par = rep(0, max(do.call('c', par_index)))

par[par_index$beta] = c(0.25, -2, 2, -0.25) 
par[par_index$alpha_tilde] = c(-1,  1, 0, 0,
                               7, -7, 0, 0,
                               -7,  7, 0, 0,
                               1, -1, 0, 0)
par[par_index$upsilon] = c(diag(c(0.25, 0.25,  4,  4,
                                  2.25, 2.25, 25, 25,
                                  2.25, 2.25, 25, 25,
                                  0.25, 0.25,  4,  4)))
par[par_index$A] = rep(0, 4)
par[par_index$R] = c(diag(c(4, 16, 16, 4)))
#    transitions:          1->2,         1->4,         2->3,         2->4, 
#                          3->1,         3->2,         3->4,         4->2, 
#                          4->5,         5->1,         5->2,         5->4
par[par_index$zeta] = c(-7.2405, 2.5, -6.2152,   1, -2.6473,  -1, -6.1475,  -1, 
                        -9.4459,  -1, -7.2404, 2.5, -7.2151,   1, -7.1778, 2.5, 
                        -5.2151,   0, -9.4459,  -1, -7.2404, 2.5, -5.2151,   0)
# par[par_index$zeta] = c(-7.2405, 2.5, -5.2152,   1, -2.6473,  -1, -5.1475,  -1, 
#                         -9.4459,  -1, -7.2404, 2.5, -5.2151,   1, -7.1778, 2.5, 
#                         -2.6523,   0, -9.4459,  -1, -7.2404, 2.5, -5.2151,   1)
par[par_index$init] = c(0, 0, 0, 0)
par[par_index$omega_tilde]= 2 * c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                                  -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                                  -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                                  -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
par[par_index$eta_omega] = rep(-1, length(par_index$eta_omega))
par[par_index$G] = c(diag(c(8, 32, 32, 8)))

A = list()
W = list()
B = list()

chosen_seed = 5
chosen_it = 21
load(paste0('Model_out/mcmc_out_1_', chosen_seed, 'it', chosen_it, '.rda'))

par = mcmc_out$chain[nrow(mcmc_out$chain), ]

b_chain = rep(1, nrow(data_format))

for(ii in 1:length(EIDs)){
    i = EIDs[ii]
    
    A[[ii]] = matrix(par[par_index$alpha_tilde], ncol =1)
    W[[ii]] = matrix(par[par_index$omega_tilde], ncol =1)
    B[[ii]] = matrix(b_chain[data_format[,"EID"] == i], ncol = 1)
}

rm(mcmc_out)

# -----------------------------------------------------------------------------

steps  = 500000
burnin =  0

s_time = Sys.time()
mcmc_out = mcmc_routine(steps, burnin, seed_num, trialNum, simulation, max_ind,
                        par, par_index, Y, x, z, B, A, W, Dn_omega, 
                        bleed_indicator)
e_time = Sys.time() - s_time; print(e_time) 