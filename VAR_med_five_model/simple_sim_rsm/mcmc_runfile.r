source('mcmc_routine.r')

args = commandArgs(TRUE)
seed_num = as.numeric(args[1])

# for(bt1 in c(1,0)) {
    bt1 = 1
    set.seed(seed_num)

    before_t1 = as.logical(bt1) # Do we handle the state changes before t1 or not?

    # Load data --------------------------------------------------------------------
    load(paste0('Data/data_format', seed_num, '.rda'))

    y = data_format[,c("y1", "y2", "y3", "y4"), drop = F]
    ids = data_format[,"id"]
    EIDs = unique(data_format[,"id"])

    # Parameter initialization ----------------------------------------------------
    par_index = list()
    par_index$alpha = 1:12
    par_index$zeta = 13:16
    par_index$diag_R = 17:20
    par_index$init = 21:22
    
    par = rep(0, max(do.call('c', par_index)))
    par[par_index$alpha] = c( 50, -3,  3,
                             100,  5, -5,
                             100, -5,  5,
                              50,  3, -3)
    par[par_index$zeta] = c(-2, -1, -1.5, -1.5)
    par[par_index$diag_R] = c(2, 2, 2, 2)
    par[par_index$init] = c(0, 0)

    n_state = 3

    B = list()
    for(i in EIDs){
        B[[i]] = matrix(data_format[data_format[,"id"] == i,"state"], ncol = 1)
    }
    # -----------------------------------------------------------------------------

    steps  = 10000
    burnin =  2000

    s_time = Sys.time()

    mcmc_out = mcmc_routine(par, par_index, B, y, ids, steps, burnin, seed_num, before_t1)

    e_time = Sys.time() - s_time; print(e_time)    
# }