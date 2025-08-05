source('mcmc_routine.r')

args = commandArgs(TRUE)
seed_num = as.numeric(args[1])
set.seed(seed_num)

dgm = T # fit the data generating model (dgm) or the approx. model

# Load data --------------------------------------------------------------------
load(paste0('Data/sim_data_', seed_num, '.rda'))

Y = data_format[, c('EID', 'hemo', 'hr', 'map', 'lactate')]
EIDs = unique(data_format[,'EID'])

# Parameter initialization -----------------------------------------------------
if(dgm) {
    par_index = list()
    par_index$alpha_tilde = 1:20
    par_index$upsilon = 21:420
    par_index$A = 421:424
    par_index$R = 425:440
    par_index$zeta = 441:452
    par_index$init = 453:456

    par = rep(0, max(do.call('c', par_index)))
    par[par_index$alpha_tilde] = c( 50, -5,  5, -2,  2,
                                   100, 10,-10,  2, -2,
                                   100,-10, 10,  2, -2,
                                    50,  5, -5, -2,  2)
    par[par_index$upsilon] = c(diag(c(25, 4, 4, 1, 1, 
                                      25, 4, 4, 1, 1, 
                                      25, 4, 4, 1, 1, 
                                      25, 4, 4, 1, 1)))
    par[par_index$A] = rep(0, 4)
    par[par_index$R] = c(diag(c(9, 9, 9, 9)))
    #    transitions:          1->2,    1->4,    2->3,    2->4, 
    #                          3->1,    3->2,    3->4,    4->2, 
    #                          4->5,    5->1,    5->2,    5->4
    par[par_index$zeta] = c(-3.7405, -4.2152, -2.6473, -2.1475, 
                            -3.4459, -2.9404, -3.2151, -3.1778, 
                            -2.0523, -3.4459, -3.2404, -3.2151)
    par[par_index$init] = c(0,0,0,0)
} else {
    par_index = list()
    par_index$alpha_tilde = 1:16
    par_index$upsilon = 17:272
    par_index$A = 273:276
    par_index$R = 277:292
    par_index$zeta = 293:304
    par_index$init = 305:308
    par_index$G = 309:324

    par = rep(0, max(do.call('c', par_index)))
    par[par_index$alpha_tilde] = c( -5,  5, -2,  2,
                                    10,-10,  2, -2,
                                   -10, 10,  2, -2,
                                     5, -5, -2,  2)
    par[par_index$upsilon] = c(diag(c(4, 4, 1, 1, 
                                      4, 4, 1, 1, 
                                      4, 4, 1, 1, 
                                      4, 4, 1, 1)))
    par[par_index$A] = rep(0, 4)
    par[par_index$R] = c(diag(c(sqrt(3), sqrt(3), sqrt(3), sqrt(3))))
    #    transitions:          1->2,    1->4,    2->3,    2->4, 
    #                          3->1,    3->2,    3->4,    4->2, 
    #                          4->5,    5->1,    5->2,    5->4
    par[par_index$zeta] = c(-3.7405, -4.2152, -2.6473, -2.1475, 
                            -3.4459, -2.9404, -3.2151, -3.1778, 
                            -2.0523, -3.4459, -3.2404, -3.2151)
    par[par_index$init] = c(0,0,0,0)
    par[par_index$G] = c(diag(c(sqrt(3), sqrt(3), sqrt(3), sqrt(3))))
}

load(paste0('Data/alpha_i_mat_', seed_num, '.rda'))
b_chain = data_format[, "b_true"]

n_state = 5

A = list()
B = list()
for(ii in 1:length(EIDs)){
    i = EIDs[ii]
    
    if(dgm) {
        A[[ii]] = alpha_i_mat[[ii]]
    } else {
        A[[ii]] = alpha_i_mat[[ii]][-c(1,6,11,16),,drop=F] # Remove baseline
    }
    
    B[[ii]] = matrix(b_chain[data_format[,"EID"] == i], ncol = 1)
    B[[ii]][1,] = 0 # we don't care about the first state
}
# -----------------------------------------------------------------------------

steps = 10000
burnin = 5000

s_time = Sys.time()

mcmc_out = mcmc_routine(steps, burnin, seed_num, par, par_index, Y, B, A, dgm)

e_time = Sys.time() - s_time; print(e_time) 