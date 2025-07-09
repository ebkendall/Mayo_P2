source('mcmc_routine.r')

args = commandArgs(TRUE)
seed_num = as.numeric(args[1])
set.seed(seed_num)

dgm = TRUE # fit the data generating model (dgm) or the approx. model

# Load data ----------------------------------------------------------------
load(paste0('Data/data_format', seed_num, '.rda'))

y = data_format[,c("y1", "y2", "y3", "y4"), drop = F]
ids = data_format[,"id"]
EIDs = unique(data_format[,"id"])

first_ind = c(0, which(diff(data_format[,"id"]) != 0)) + 1

# Parameter initialization -------------------------------------------------
if(dgm) {
    par_index = list()
    par_index$alpha = 1:12
    par_index$zeta = 13:16
    par_index$diag_R = 17:20
    par_index$init = 21:22
    
    par = rep(0, max(do.call('c', par_index)))
    par[par_index$alpha] = c( 50,  -5,   5,
                             100,  10, -10,
                             100, -10,  10,
                              50,   5,  -5)
    par[par_index$zeta] = c(-2, -2, -1.5, -1.5)
    par[par_index$diag_R] = c(1.386294, 1.386294, 1.386294, 1.386294)
    par[par_index$init] = c(0, 0)
} else {
    par_index = list()
    par_index$alpha = 1:8
    par_index$zeta = 9:12
    par_index$diag_R = 13:16
    par_index$init = 17:18
    par_index$diag_G = 19:22
    
    par = rep(0, max(do.call('c', par_index)))
    par[par_index$alpha] = c( -5,   5,
                              10, -10,
                             -10,  10,
                               5,  -5)
    par[par_index$zeta] = c(-2, -2, -1.5, -1.5)
    par[par_index$diag_R] = c(1.386294, 1.386294, 1.386294, 1.386294)
    par[par_index$init] = c(0, 0)
    par[par_index$diag_G] = c(1.386294, 1.386294, 1.386294, 1.386294)
}

n_state = 3

B = list()
for(i in EIDs){
    B[[i]] = matrix(data_format[data_format[,"id"] == i,"state"], ncol = 1)
    B[[i]][1,] = 0 # we don't care about the first state
}
# -----------------------------------------------------------------------------

steps  = 10000
burnin =  2000

s_time = Sys.time()

mcmc_out = mcmc_routine(par, par_index, B, y, ids, steps, burnin, seed_num, dgm)

e_time = Sys.time() - s_time; print(e_time)    