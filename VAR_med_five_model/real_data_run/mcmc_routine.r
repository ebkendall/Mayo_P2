library(mvtnorm, quietly=T)
library(Matrix, quietly=T)
library(LaplacesDemon, quietly=T)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(expm, quietly = T)
sourceCpp("likelihood_fnc.cpp")

# Needed for OpenMP C++ parallel
# Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
# Sys.setenv("PKG_LIBS" = "-fopenmp")

mcmc_routine = function( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, 
                         trialNum, Dn_omega, simulation, bleed_indicator, max_ind, 
                         df_num, sampling_num, states_per_step, steps_per_it){
    
    EIDs = as.character(unique(Y[,'EID']))
    
    # Number of cores over which to parallelize --------------------------------
    n_cores = 10#strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
    print(paste0("Number of cores: ", n_cores))
    
    # Transition information ---------------------------------------------------
    adj_mat = matrix(c(1, 1, 0, 1, 0,
                       0, 1, 1, 1, 0,
                       1, 1, 1, 1, 0,
                       0, 1, 0, 1, 1,
                       1, 1, 0, 1, 1), nrow=5, byrow = T)
    adj_mat_sub = matrix(c(1, 0, 0, 1, 0,
                           0, 1, 0, 0, 0,
                           1, 0, 1, 1, 0,
                           0, 0, 0, 1, 1,
                           1, 0, 0, 1, 1), nrow=5, byrow = T)
    initialize_cpp(adj_mat, adj_mat_sub, states_per_step)

    # Index of observed versus missing data ------------------------------------
    # 1 = observed, 0 = missing
    otype = !is.na(Y[, c('hemo','hr','map','lactate')])
    colnames(otype) = c('hemo','hr','map','lactate')

    # Metropolis Parameter Index for MH within Gibbs updates -------------------
    mpi = list(c(par_index$vec_init),
               c(par_index$vec_zeta[c(1,5, 7,11,15,21)]), # baselines (w/    S2)
               c(par_index$vec_zeta[c(3,9,13,17,19,23)]), # baselines (w/out S2)
               c(par_index$vec_zeta[c(2,12,16,22)]),      # RBC > 0 (to S2)
               c(par_index$vec_zeta[c(4,14,18,24)]),      # RBC > 0 (no S2)
               c(par_index$vec_zeta[c(6,8,10,20)]),       # RBC < 0
               c(par_index$vec_A[c(1,5,9,13,17)]),
               c(par_index$vec_A[c(2,6,10,14,18)]),
               c(par_index$vec_A[c(3,7,11,15,19)]),
               c(par_index$vec_A[c(4,8,12,16,20)]),
               c(par_index$vec_upsilon_omega[c(1:16, 36:57)]),
               c(par_index$vec_upsilon_omega[c(17:35, 58:88)]),
               c(par_index$vec_R))

    n_group = length(mpi)

    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))
    pscale = rep( 0.0001, n_group)
  
    # Initialize data storage --------------------------------------------------
    chain_length_MASTER = 10000
    chain = matrix( 0, chain_length_MASTER, length(par)) 
    B_chain = matrix( 0, chain_length_MASTER, nrow(Y)) 
    hc_chain = hr_chain = bp_chain = la_chain = NULL
    accept = rep( 0, n_group)
    
    # Initialize design matrix -------------------------------------------------
    Xn = initialize_X( as.numeric(EIDs), Y, x)
    
    # Initialize Y and states --------------------------------------------------
    if(!simulation) {
        vital_means = colMeans(Y[,c('hemo', 'hr', 'map', 'lactate')], na.rm = T)
        Y_B_Dn_init = initialize_Y(as.numeric(EIDs), par, par_index, A, Y, z, Xn,
                                   Dn_omega, W, otype, n_cores, vital_means)
        
        Y[, c('hemo', 'hr', 'map', 'lactate')] = Y_B_Dn_init[[1]]
        B = Y_B_Dn_init[[2]]; names(B) = EIDs
        Dn = Y_B_Dn_init[[3]]; names(Dn) = EIDs
        
        for(i in 1:100) {
            Y = update_Y_i_cpp( as.numeric(EIDs), par, par_index, A, Y, Dn, Xn,
                                otype, Dn_omega, W, B, n_cores)
            colnames(Y) = c('EID','hemo', 'hr', 'map', 'lactate',
                            'RBC_rule', 'clinic_rule')    
        }    
    } else {
        B_Dn = mle_state_seq(as.numeric(EIDs), par, par_index, A, Y, z, Xn, 
                             Dn_omega, W, n_cores)
        B = B_Dn[[1]]; names(B) = EIDs
        Dn = B_Dn[[2]]; names(Dn) = EIDs    
    }
    
    # Keeping track of the sampled alpha_i -------------------------------------
    A_chain = vector(mode = "list", length = 5)
    a_chain_id = c(3, 86, 163, 237, 427)
    for(a_ind in 1:length(A_chain)) {
        A_chain[[a_ind]] = matrix(nrow = length(par_index$vec_alpha_tilde), 
                                  ncol = chain_length_MASTER)
    }

    # Start Metropolis-within-Gibbs Algorithm ----------------------------------
    chain[1,] = par
    B_chain[1, ] = do.call( 'c', B)
    
    mcmc_start_t = Sys.time()
    for(ttt in 2:steps){

        # Every 10,000 steps, save existing MCMC results -----------------------
        if(ttt %% chain_length_MASTER == 0) {
            chain_ind = chain_length_MASTER
        } else {
            chain_ind = ttt - chain_length_MASTER * floor(ttt / chain_length_MASTER)
        }

        # Imputing the missing Y values ----------------------------------------
        if(!simulation) {
            Y = update_Y_i_cpp( as.numeric(EIDs), par, par_index, A, Y, Dn, Xn,
                                otype, Dn_omega, W, B, n_cores)
            colnames(Y) = c('EID','hemo', 'hr', 'map', 'lactate',
                            'RBC_rule', 'clinic_rule')
        }

        bbb_start_t = Sys.time()
        for(s in 1:steps_per_it) {
            # Random sample update ---------------------------------------------
            if(sampling_num == 1) {
                B_Dn = mh_up(as.numeric(EIDs), par, par_index, A, B, Y, 
                             z, Xn, Dn_omega, W, bleed_indicator, 
                             n_cores, states_per_step);
                B = B_Dn[[1]]; names(B) = EIDs
                Dn = B_Dn[[2]]; names(Dn) = EIDs
            }
            
            # Almost-Gibbs update ----------------------------------------------
            else if(sampling_num == 2) {
                B_Dn = almost_gibbs_up(as.numeric(EIDs), par, par_index, A, B, Y, 
                                       z, Xn, Dn_omega, W, bleed_indicator, 
                                       n_cores, states_per_step);
                B = B_Dn[[1]]; names(B) = EIDs
                Dn = B_Dn[[2]]; names(Dn) = EIDs
            } 
            
            # Gibbs update -----------------------------------------------------
            else if(sampling_num == 3) {
                B_Dn = gibbs_up(as.numeric(EIDs), par, par_index, A, B, Y, z, Xn, 
                                Dn_omega, W, bleed_indicator, n_cores, states_per_step);
                B = B_Dn[[1]]; names(B) = EIDs
                Dn = B_Dn[[2]]; names(Dn) = EIDs
            }
            
            # Full seq MH update -----------------------------------------------
            else if(sampling_num == 4) {
                B_Dn = mh_up_all(as.numeric(EIDs), par, par_index, A, B, Y, z, Xn, 
                                 Dn_omega, W, bleed_indicator, n_cores);
                B = B_Dn[[1]]; names(B) = EIDs
                Dn = B_Dn[[2]]; names(Dn) = EIDs
            }
        }
        bbb_end_t = Sys.time() - bbb_start_t; print(bbb_end_t)
        
        # Gibbs: alpha_i -------------------------------------------------------
        A = update_alpha_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn,
                                Dn_omega, W, B, n_cores)
        names(A) = EIDs

        for(aaa in 1:length(a_chain_id)) {
            A_chain[[aaa]][,chain_ind] = A[[a_chain_id[aaa]]]
        }

        # Gibbs: omega_i -------------------------------------------------------
        W = update_omega_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn,
                                Dn_omega, A, B, n_cores)
        names(W) = EIDs

        # # Gibbs: alpha~, omega~, beta, & Upsilon -------------------------------
        # par = update_alpha_tilde_cpp( as.numeric(EIDs), par, par_index, A, Y)
        # par = update_omega_tilde_cpp( as.numeric(EIDs), par, par_index, W, Y)
        # par = update_beta_upsilon_cpp( as.numeric(EIDs), par, par_index, A, Y,
        #                                  Dn, Xn, Dn_omega, W, B, n_cores)

        # Store current parameter updates --------------------------------------
        chain[chain_ind,] = par
        
        # # Evaluate log-likelihood before MH step -------------------------------
        # log_target_prev = log_post_cpp( as.numeric(EIDs), par, par_index, A, B,
        #                                 Y, z, Dn, Xn, Dn_omega, W, n_cores)
        # 
        # if(!is.finite(log_target_prev)){
        #     print("Infinite log-posterior; Gibbs update went wrong")
        #     print(paste0("value: ", log_target_prev))
        #     break
        # }
        # 
        # # Metropolis-Hastings updates ------------------------------------------
        # for(j in 1:n_group) {
        # 
        #     ind_j = mpi[[j]]
        #     proposal = par
        # 
        #     # MH update: A, zeta, init -----------------------------------------
        #     if(sum(ind_j %in% par_index$vec_R) == 0) {
        #         proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
        #                                    sigma=pscale[[j]]*pcov[[j]])
        # 
        #         # Evaluate proposed log-likelihood -----------------------------
        #         log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index,
        #                                    A, B, Y, z, Dn, Xn, Dn_omega, W, n_cores)
        # 
        #         if(ttt < burnin){
        #             while(!is.finite(log_target)){
        #                 print('bad proposal')
        #                 proposal = par
        #                 proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
        #                                            sigma=pcov[[j]]*pscale[j])
        # 
        #                 log_target = log_post_cpp( as.numeric(EIDs), proposal,
        #                                            par_index, A, B, Y, z, Dn,
        #                                            Xn, Dn_omega, W, n_cores)
        #             }
        #         }
        # 
        #         if(!is.finite(log_target) | is.nan(log_target)) {
        #             # Ensuring that we do not have problems from C++
        #             print(paste0("bad proposal post burnin: ", log_target))
        #             log_target = -Inf
        #         }
        # 
        #         # Compute Metropolis ratio -------------------------------------
        #         if( log_target - log_target_prev > log(runif(1,0,1)) ){
        #             log_target_prev = log_target
        #             par[ind_j] = proposal[ind_j]
        #             accept[j] = accept[j] +1
        #         }
        # 
        #         chain[chain_ind,ind_j] = par[ind_j]
        # 
        #         # Proposal tuning scheme ---------------------------------------
        #         # During the burnin period, update the proposal covariance in each step
        #         # to capture the relationships within the parameters vectors for each
        #         # transition.  This helps with mixing.
        #         if(ttt < burnin){
        #             if(ttt == 100)  pscale[j] = 1
        #             
        #             if(100 <= ttt & ttt <= 2000){
        #                 temp_chain = chain[1:ttt,ind_j]
        #                 pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
        #                 
        #             } else if(2000 < ttt){
        #                 temp_chain = chain[(ttt-2000):ttt,ind_j]
        #                 pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
        #             }
        #             if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )
        #             
        #             # Tune the proposal covariance for each transition to achieve
        #             # reasonable acceptance ratios.
        #             if(ttt %% 30 == 0){
        #                 if(ttt %% 480 == 0){
        #                     accept[j] = 0
        #                     
        #                 } else if( accept[j] / (ttt %% 480) < .4 ){ 
        #                     pscale[j] = (.75^2)*pscale[j]
        #                     
        #                 } else if( accept[j] / (ttt %% 480) > .5 ){ 
        #                     pscale[j] = (1.25^2)*pscale[j]
        #                 }
        #             }
        #         }
        #         # --------------------------------------------------------------
        # 
        #     } else {
        #         # MH update: R -------------------------------------------------
        #         # Prior for R 
        #         nu_R = 1000
        #         psi_R = diag(c(9, 9, 9, 9))
        #         psi_R = (nu_R - 4 - 1) * psi_R
        #         
        #         # Proposal df and scale matrix
        #         psi_nu_q = proposal_R_cpp_new(nu_R, psi_R, Y, Dn, Xn,
        #                                       A, par, par_index, as.numeric(EIDs),
        #                                       B, Dn_omega, W)
        #         psi_q = psi_nu_q[[1]]
        #         nu_q = psi_nu_q[[2]]
        #         
        #         # Current
        #         curr_R = matrix(par[ind_j], nrow = 4)
        #         log_q_curr = dinvwishart(Sigma = curr_R, nu = nu_q, S = psi_q, log = T)
        #         
        #         # Proposal
        #         proposal[ind_j] = c(rinvwishart(nu = nu_q, S = psi_q))
        #         prop_R = matrix(proposal[ind_j], nrow = 4)
        #         log_q_prop = dinvwishart(Sigma = prop_R, nu = nu_q, S = psi_q, log = T)
        #         
        #         log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index,
        #                                    A, B, Y, z, Dn, Xn, Dn_omega, W, n_cores)
        #         
        #         if(!is.finite(log_target + log_q_curr - log_target_prev - log_q_prop) | 
        #            is.nan(log_target + log_q_curr - log_target_prev - log_q_prop)) {
        #             print("Current R")
        #             print(curr_R)
        #             print("Prop R")
        #             print(prop_R)
        #             print(paste0("log_target: ", log_target))
        #             print(paste0("log_q_curr: ", log_q_curr))
        #             print(paste0("log_target_prev: ", log_target_prev))
        #             print(paste0("log_q_prop: ", log_q_prop))
        #             print(log_target + log_q_curr - log_target_prev - log_q_prop)
        #         }
        #         
        #         # Compute Metropolis ratio 
        #         if( log_target + log_q_curr - log_target_prev - log_q_prop > log(runif(1,0,1)) ){
        #             log_target_prev = log_target
        #             par[ind_j] = proposal[ind_j]
        #             accept[j] = accept[j] +1
        #         }
        #         
        #         chain[chain_ind,ind_j] = par[ind_j]
        #     }
        # }

        # Restart the acceptance ratio at burnin
        if(ttt == burnin) accept = rep( 0, n_group)
        
        B_chain[ chain_ind, ] = do.call( 'c', B)
        # hc_chain[chain_ind, ] = Y[,'hemo']
        # hr_chain[chain_ind, ] = Y[,'hr']
        # bp_chain[chain_ind, ] = Y[,'map']
        # la_chain[chain_ind, ] = Y[,'lactate']
        # ----------------------------------------------------------------------

        if(ttt%%1==0)  cat('--->',ttt,'\n')
        # if(ttt%%1==0) print(log_target_prev)
        if(ttt%%100==0) print(accept)
        
        if(ttt > burnin & ttt%%chain_length_MASTER==0) {
            mcmc_end_t = Sys.time() - mcmc_start_t; print(mcmc_end_t)
            index_keep = seq(1, chain_length_MASTER, by = 5)
            mcmc_out_temp = list(chain    = chain[index_keep,], 
                                 B_chain  = B_chain[index_keep,], 
                                 hc_chain = Y[,'hemo'],#hc_chain[index_keep,], 
                                 hr_chain = Y[,'hr'], #hr_chain[index_keep,],
                                 bp_chain = Y[,'map'], #bp_chain[index_keep,], 
                                 la_chain = Y[,'lactate'], #la_chain[index_keep,],
                                 A_chain  = A_chain,
                                 otype=otype, accept=accept/length(burnin:ttt), 
                                 pscale=pscale, pcov = pcov, par_index=par_index)
            if(simulation) {
                save(mcmc_out_temp, file = paste0('Model_out/mcmc_out_',ind,'_', 
                                                  trialNum, 'it', ttt/chain_length_MASTER + (max_ind - 5), 
                                                  '_samp', sampling_num, '_', states_per_step, '_', steps_per_it,
                                                  '_sim.rda'))
            } else {
                save(mcmc_out_temp, file = paste0('Model_out/mcmc_out_',ind,'_', 
                                                  trialNum, 'it', ttt/chain_length_MASTER + (max_ind - 5), 
                                                  '.rda'))
            }
            # Reset the chains
            chain = matrix( NA, chain_length_MASTER, length(par)) 
            B_chain = hc_chain = hr_chain = bp_chain = la_chain = matrix( NA, chain_length_MASTER, nrow(Y))
        }
    }
    # ---------------------------------------------------------------------------

    return(list( chain=chain, B_chain=B_chain, hc_chain=hc_chain, hr_chain=hr_chain, 
                bp_chain=bp_chain, la_chain = la_chain, A_chain = A_chain, otype=otype,
                accept=accept/(steps-burnin), pscale=pscale, pcov = pcov,
                par_index=par_index))
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
            e_time_24 = s_time + 1440
            RBC_diff_12 = RBC_diff_24 = 0
            
            if (e_time_12 <= max_time) {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                ind_12 = order(abs(sub_dat[,"time"] - e_time_12))[1]
                RBC_diff_12 = sub_dat[ind_12, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            } else {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
                RBC_diff_12 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            }
            if (e_time_24 <= max_time) {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                ind_24 = order(abs(sub_dat[,"time"] - e_time_24))[1]
                RBC_diff_24 = sub_dat[ind_24, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            } else {
                s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
                e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
                RBC_diff_24 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
            }
            
            if(RBC_diff_12 >=3 | RBC_diff_24 >= 6) {
                admin_times = sub_dat[sub_dat[,"RBC_admin"] != 0, "time"]
                if(RBC_diff_12 >=3) {
                    a_t = which(admin_times >= s_time & admin_times < e_time_12)
                    first_time = admin_times[a_t[1]]
                    order_times = sub_dat[sub_dat[,"RBC_ordered"] != 0, "time"]
                    if(sum(order_times <= first_time) == 0) {
                        # print(paste0(i, ", ", sub_dat[1,"EID"]))
                        first_order_time = first_time
                    } else {
                        first_order_time = max(order_times[order_times <= first_time])   
                    }
                } else if (RBC_diff_24 >= 6) {
                    a_t = which(admin_times >= s_time & admin_times < e_time_24)
                    first_time = admin_times[a_t[1]]  
                    order_times = sub_dat[sub_dat[,"RBC_ordered"] != 0, "time"]
                    if(sum(order_times <= first_time) == 0) {
                        # print(paste0(i, ", ", sub_dat[1,"EID"]))
                        first_order_time = first_time
                    } else {
                        first_order_time = max(order_times[order_times <= first_time])   
                    }
                }
                
                bleed_indicator[data_format[,"EID"] == bleed_pat[i] & 
                                    data_format[,"time"] == first_order_time] = 1
                break
            }
            
        }
    }
    return(bleed_indicator)
}


var_R_calc <- function(psi, nu, p) {
    var_mat = matrix(0, p, p)
    for(i in 1:p) {
        for(j in 1:p) {
            num_ij = (nu - p + 1) * (psi[i,j] * psi[i,j]) + (nu - p - 1) * psi[i,i] * psi[j,j]
            den_ij = (nu - p) * (nu - p - 1) * (nu - p - 1) * (nu - p - 3)
            var_mat[i,j] = num_ij / den_ij
        }
    }
    return(var_mat)
}