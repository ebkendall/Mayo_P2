source('test_set_routine.r')

args = commandArgs(TRUE)
seed_num = as.numeric(args[1])
set.seed(seed_num)

simulation = F

# Load data --------------------------------------------------------------------
load('Data/data_format_TEST.rda')

# Remove clinic and RBC rule to leave this as uninformed
data_format[,'RBC_rule'] = 0
data_format[,'clinic_rule'] = 0
bleed_indicator = rep(0, nrow(data_format))

load('Data/Dn_omega_TEST.rda')

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

load('Model_out/post_med_stack_it4_0.rda')
par = post_med_stack

A = list()
B = list()

b_chain = rep(1, nrow(data_format))

for(ii in 1:length(EIDs)){
    i = EIDs[ii]
    
    A[[ii]] = matrix(par[par_index$alpha_tilde], ncol =1)
    B[[ii]] = matrix(b_chain[data_format[,"EID"] == i], ncol = 1)
}

# -----------------------------------------------------------------------------

steps  = 500000
burnin =  0

s_time = Sys.time()
mcmc_out = mcmc_routine(steps, burnin, seed_num, trialNum, simulation, max_ind,
                        par, par_index, Y, x, z, B, A, Dn_omega, bleed_indicator)
e_time = Sys.time() - s_time; print(e_time) 