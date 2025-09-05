library(mvtnorm, quietly=T)
library(LaplacesDemon, quietly=T)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("likelihood_fnc.cpp")

# Needed for OpenMP C++ parallel
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
Sys.setenv("PKG_LIBS" = "-fopenmp")

mcmc_routine = function(steps, burnin, seed_num, trialNum, simulation, max_ind,
                        par, par_index, Y, x, z, B, A, W, Dn_omega, 
                        bleed_indicator) {
    
    EIDs = as.numeric(unique(Y[,'EID']))
    
    # Number of cores over which to parallelize --------------------------------
    n_cores = 1
    print(paste0("Number of cores: ", n_cores))
    
    # Transition information ---------------------------------------------------
    adj_mat = matrix(c(1, 1, 0, 1, 0,
                       0, 1, 1, 1, 0,
                       1, 1, 1, 1, 0,
                       0, 1, 0, 1, 1,
                       1, 1, 0, 1, 1), nrow=5, byrow = T)
    adj_mat_sub = matrix(c(1, 0, 0, 1, 0,
                           0, 1, 0, 0, 0,
                           0, 0, 1, 0, 0,
                           0, 0, 0, 1, 1,
                           1, 0, 0, 1, 1), nrow=5, byrow = T)
    initialize_cpp(adj_mat, adj_mat_sub)
    
    # Index of observed versus missing data ------------------------------------
    # 1 = observed, 0 = missing
    otype = !is.na(Y[, c('hemo','hr','map','lactate')])
    colnames(otype) = c('hemo','hr','map','lactate')
    
    # Initialize design matrix -------------------------------------------------
    Xn = initialize_Xn(EIDs, Y, x)
    
    # Initialize Y and states --------------------------------------------------
    impute_step = TRUE
    if(sum(c(otype) == 0) > 0) {
        
        # There exists missing Y values
        print("Missing Y values - initializing")
        
        vital_means = colMeans(Y[,c('hemo', 'hr', 'map', 'lactate')], na.rm = T)
        Y_init = initialize_Y(EIDs, par, par_index, A, W, Y, z, Dn_omega, 
                              Xn, otype, n_cores, vital_means)
        
        Y[, c('hemo', 'hr', 'map', 'lactate')] = Y_init
        
        B_Dn = mle_state_seq(EIDs, par, par_index, A, W, Y, z, 
                             Dn_omega, Xn, n_cores)
        B = B_Dn[[1]]
        Dn_alpha = B_Dn[[2]]
        
        print("Done")
        
    } else {
        
        # There are NO missing Y values
        print("No missingness")
        impute_step = FALSE
        
        B_Dn = mle_state_seq(EIDs, par, par_index, A, W, Y, z, 
                             Dn_omega, Xn, n_cores)
        B = B_Dn[[1]]
        Dn_alpha = B_Dn[[2]]
        
    }
    
    # Initialize data storage --------------------------------------------------
    reset_step = 10000
    chain_length_MASTER = 1000
    
    chain = matrix(NA, reset_step, length(par)) 
    B_chain = matrix(NA, chain_length_MASTER, nrow(Y)) 
    
    # Start Metropolis-within-Gibbs Algorithm ----------------------------------
    chain[1,] = par
    B_chain[1, ] = do.call( 'c', B)
    
    if(impute_step) {
        hc_chain = hr_chain = matrix(NA, chain_length_MASTER, nrow(Y))
        bp_chain = la_chain = matrix(NA, chain_length_MASTER, nrow(Y))   
        
        hc_chain[1, ] = Y[,'hemo']
        hr_chain[1, ] = Y[,'hr']
        bp_chain[1, ] = Y[,'map']
        la_chain[1, ] = Y[,'lactate']    
    }
    
    mcmc_start_t = Sys.time()
    for(ttt in 2:steps){
        ttt_start_t = Sys.time()
        
        # Every 10,000 steps, save existing MCMC results -----------------------
        chain_ind = NULL
        chain_ttt = NULL
        if(ttt %% reset_step == 0) {
            chain_ind = chain_length_MASTER
            chain_ttt = reset_step
        } else {
            new_t = ttt - floor(ttt / reset_step) * reset_step
            chain_ind = floor((new_t - 1)/10) + 1
            chain_ttt = ttt %% reset_step
        }
        
        # Sample gamma_1 -------------------------------------------------------
        gamma_1 = update_gamma_i(EIDs, par, par_index, A, W, Y, 
                                 Dn_alpha, Dn_omega, Xn, n_cores)
        
        # Imputing the missing Y values ----------------------------------------
        if(impute_step) {
            Y = impute_Y(EIDs, par, par_index, A, W, Y, Dn_alpha, Dn_omega, Xn, 
                         gamma_1, otype, n_cores)
            colnames(Y) = c('EID','hemo', 'hr', 'map', 'lactate',
                            'RBC_rule', 'clinic_rule')    
        }
        
        # State-space update (B) -----------------------------------------------
        sps = sample(x = 2:50, size = 1, replace = T) # sps >= 2
        B_Dn = state_sampler(EIDs, par, par_index, B, A, W, Y, z, Dn_alpha, 
                             Dn_omega, Xn, bleed_indicator, gamma_1, sps, n_cores)
        B = B_Dn[[1]]
        Dn_alpha = B_Dn[[2]]
        
        # Gibbs: alpha_i -------------------------------------------------------
        A = update_alpha_i(EIDs, par, par_index, W, Y, Dn_alpha, 
                           Dn_omega, Xn, gamma_1, n_cores)
        
        # Gibbs: omega_i -------------------------------------------------------
        W = update_omega_i(EIDs, par, par_index, A, Y, Dn_alpha, 
                           Dn_omega, Xn, gamma_1, n_cores)
        
        # Store information ----------------------------------------------------
        B_chain[ chain_ind, ] = do.call( 'c', B)
        
        if(impute_step) {
            hc_chain[chain_ind, ] = Y[,'hemo']
            hr_chain[chain_ind, ] = Y[,'hr']
            bp_chain[chain_ind, ] = Y[,'map']
            la_chain[chain_ind, ] = Y[,'lactate']
        }
        # ----------------------------------------------------------------------
        
        chain[chain_ttt,] = par
        
        cat('--> ', ttt, ', c_ttt = ', chain_ttt, ', c_ind = ', chain_ind, '\n')
        
        ttt_end_t = Sys.time() - ttt_start_t; print(ttt_end_t)
        
        if(ttt > burnin & ttt%%reset_step==0) {
            mcmc_end_t = Sys.time() - mcmc_start_t; print(mcmc_end_t)
            
            index_keep = seq(10, reset_step, by = 10)
            index_keep_2 = seq(2, chain_length_MASTER, by = 2)
            
            mcmc_out = list(chain    = chain[index_keep,], 
                            B_chain  = B_chain[index_keep_2,], 
                            hc_chain = hc_chain[index_keep_2,], 
                            hr_chain = hr_chain[index_keep_2,],
                            bp_chain = bp_chain[index_keep_2,], 
                            la_chain = la_chain[index_keep_2,],
                            alpha_i = A,
                            omega_i = W,
                            otype=otype, par_index=par_index)
            
            save(mcmc_out, 
                 file = paste0('Model_out/mcmc_out_1_', seed_num, 
                               'it', ttt/reset_step, '_TEST.rda'))
            
            # Reset the chains
            chain = matrix(NA, reset_step, length(par)) 
            B_chain = matrix(NA, chain_length_MASTER, nrow(Y)) 
            if(impute_step) {
                hc_chain = hr_chain = matrix(NA, chain_length_MASTER, nrow(Y))
                bp_chain = la_chain = matrix(NA, chain_length_MASTER, nrow(Y))    
            }
        }
    }
    # --------------------------------------------------------------------------
    
    return(0)
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

b_ind_fnc <- function(data_format) {
    bleed_pat = unique(data_format[data_format[,"RBC_rule"] != 0, "EID"])
    bleed_indicator = rep(0, nrow(data_format))
    for(i in 1:length(bleed_pat)) {
        sub_dat = data_format[data_format[,"EID"] == bleed_pat[i], ]
        
        # Check in any 12 hour period
        max_time = tail(sub_dat[,"time"], 1)
        when_rbc = c(1, which(diff(sub_dat[,"n_RBC_admin"]) != 0))
        
        for(j in 1:length(when_rbc)) {
            s_time = sub_dat[when_rbc[j], "time"]
            e_time_12 = s_time + 720
            RBC_diff_12 = 0
            
            if (e_time_12 <= max_time) {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                ind_12 = order(abs(sub_dat[,"time"] - e_time_12))[1]
                RBC_diff_12 = sub_dat[ind_12, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            } else {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
                RBC_diff_12 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            }
            
            if(RBC_diff_12 >=3) {
                admin_times = sub_dat[sub_dat[,"RBC_admin"] != 0, "time"]
                
                a_t = which(admin_times >= s_time & admin_times < e_time_12)
                
                # first_time = admin_times[a_t[1]]       # before first of three RBCs
                first_time = admin_times[tail(a_t, 1)] # before last of three RBCs
                
                order_times = sub_dat[sub_dat[,"RBC_ordered"] != 0, "time"]
                if(sum(order_times <= first_time) == 0) {
                    first_order_time = first_time
                } else {
                    first_order_time = max(order_times[order_times <= first_time])   
                }
                
                bleed_indicator[data_format[,"EID"] == bleed_pat[i] & 
                                    data_format[,"time"] == first_order_time] = 1
                break
            }
            
        }
    }
    return(bleed_indicator)
}
