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
                         df_num, sampling_num){
    
    n_cores = strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
    print(paste0("Number of cores: ", n_cores))

    EIDs = as.character(unique(Y[,'EID']))
    
    # Number of states sampled at a time ---------------------------------------
    # *********** DONT FORGET TO CHANGE THIS NUMBER IN THE .cpp FILE ***********
    t_pt_length = 3 

    # Index of observed versus missing data ------------------------------------
    # 1 = observed, 0 = missing
    otype = !is.na(Y[, c('hemo','hr','map','lactate')])
    colnames(otype) = c('hemo','hr','map','lactate')

    # Metropolis Parameter Index for MH within Gibbs updates -------------------
    mpi = list(# c(par_index$vec_init),
               c(par_index$vec_zeta[c(1,5, 7,11,15,21)]), # baselines (w/    S2)
               c(par_index$vec_zeta[c(3,9,13,17,19,23)]), # baselines (w/out S2)
               c(par_index$vec_zeta[c(2,12,16,22)]),      # RBC > 0 (to S2)
               c(par_index$vec_zeta[c(4,14,18,24)]),      # RBC > 0 (no S2)
               c(par_index$vec_zeta[c(6,8,10,20)]),       # RBC < 0 
               c(par_index$vec_A[c(1,5,9,13,17)]),
               c(par_index$vec_A[c(2,6,10,14,18)]),
               c(par_index$vec_A[c(3,7,11,15,19)]),
               c(par_index$vec_A[c(4,8,12,16,20)]),
               # c(par_index$vec_upsilon_omega[c(1:16, 36:57)]),
               # c(par_index$vec_upsilon_omega[c(17:35, 58:88)]),
               c(par_index$vec_R))

    n_group = length(mpi)

    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))*.001
    pscale = rep( 1, n_group)

    # if(!simulation) {
    # } 
  
    # Initialize data storage --------------------------------------------------
    chain_length_MASTER = 10000
    chain = matrix( 0, chain_length_MASTER, length(par)) 
    B_chain = hc_chain = hr_chain = bp_chain = la_chain = matrix( 0, chain_length_MASTER, nrow(Y)) 
    accept = rep( 0, n_group)

    # Initialize design matrices -----------------------------------------------
    Dn_Xn = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x, n_cores)
    Dn = Dn_Xn[[1]]; names(Dn) = EIDs
    Xn = Dn_Xn[[2]]
    
    # Keeping track of the sampled alpha_i -------------------------------------
    A_chain = vector(mode = "list", length = 5)
    a_chain_id = c(3, 86, 163, 237, 427)
    for(a_ind in 1:length(A_chain)) {
        A_chain[[a_ind]] = matrix(nrow = length(par_index$vec_alpha_tilde), 
                                  ncol = chain_length_MASTER)
    }

    # Start Metropolis-within-Gibbs Algorithm ----------------------------------
    chain[1,] = par
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

        # Gibbs: alpha_i -------------------------------------------------------
        A = update_alpha_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn,
                                Dn_omega, W, B, n_cores)
        names(A) = EIDs

        for(aaa in 1:length(a_chain_id)) {
            A_chain[[aaa]][,chain_ind] = A[[a_chain_id[aaa]]]
        }

        # Gibbs: omega_i -------------------------------------------------------
        # W = update_omega_i_cpp( as.numeric(EIDs), par, par_index, Y, Dn, Xn,
        #                         Dn_omega, A, B, n_cores)
        # names(W) = EIDs
        
        # Metropolis-within-Gibbs: B (states) ----------------------------------
        if(sampling_num <= 3) {
            B_Dn = update_b_i_MH(as.numeric(EIDs), par, par_index, A, B, Y, z, Dn,
                                Xn, Dn_omega, W, bleed_indicator, n_cores, 
                                t_pt_length, sampling_num)
        } else {
            B_Dn = update_b_i_gibbs(as.numeric(EIDs), par, par_index, A, B, Y, z, Dn,
                                    Xn, Dn_omega, W, bleed_indicator, n_cores, 
                                    t_pt_length)
        }
        B = B_Dn[[1]]; names(B) = EIDs
        Dn = B_Dn[[2]]; names(Dn) = EIDs
        
        # Gibbs: alpha~, omega~, beta, & Upsilon -------------------------------
        par = update_alpha_tilde_cpp( as.numeric(EIDs), par, par_index, A, Y)
        # par = update_omega_tilde_cpp( as.numeric(EIDs), par, par_index, W, Y)
        # par = update_beta_Upsilon_R_cpp( as.numeric(EIDs), par, par_index, A, Y,
        #                                  Dn, Xn, Dn_omega, W, B, n_cores)


        # Store current parameter updates --------------------------------------
        chain[chain_ind,] = par
        
        # Evaluate log-likelihood before MH step -------------------------------
        log_target_prev = log_post_cpp( as.numeric(EIDs), par, par_index, A, B,
                                        Y, z, Dn, Xn, Dn_omega, W, n_cores)
        
        if(!is.finite(log_target_prev)){
            print("Infinite log-posterior; Gibbs update went wrong")
            print(paste0("value: ", log_target_prev))
            break
        }
        
        # Metropolis-Hastings updates ------------------------------------------
        for(j in 1:n_group) {

            ind_j = mpi[[j]]
            proposal = par
        
            # MH update: A, zeta, init -----------------------------------------
            if(sum(ind_j %in% par_index$vec_R) == 0) {
                proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                           sigma=pscale[[j]]*pcov[[j]])
        
                # Evaluate proposed log-likelihood -----------------------------
                log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index,
                                           A, B, Y, z, Dn, Xn, Dn_omega, W, n_cores)
        
                if(ttt < burnin){
                    while(!is.finite(log_target)){
                        print('bad proposal')
                        proposal = par
                        proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                                   sigma=pcov[[j]]*pscale[j])
        
                        log_target = log_post_cpp( as.numeric(EIDs), proposal,
                                                   par_index, A, B, Y, z, Dn,
                                                   Xn, Dn_omega, W, n_cores)
                    }
                }
        
                if(!is.finite(log_target) | is.nan(log_target)) {
                    # Ensuring that we do not have problems from C++
                    print(paste0("bad proposal post burnin: ", log_target))
                    log_target = -Inf
                }
        
                # Compute Metropolis ratio -------------------------------------
                if( log_target - log_target_prev > log(runif(1,0,1)) ){
                    log_target_prev = log_target
                    par[ind_j] = proposal[ind_j]
                    accept[j] = accept[j] +1
                }
        
                chain[chain_ind,ind_j] = par[ind_j]
        
                # Proposal tuning scheme ---------------------------------------
                # During the burnin period, update the proposal covariance in each step
                # to capture the relationships within the parameters vectors for each
                # transition.  This helps with mixing.
                if(ttt < burnin){
                    if(ttt == 100)  pscale[j] = 1
                    
                    if(100 <= ttt & ttt <= 2000){
                        temp_chain = chain[1:ttt,ind_j]
                        pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                        
                    } else if(2000 < ttt){
                        temp_chain = chain[(ttt-2000):ttt,ind_j]
                        pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                    }
                    if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )
                    
                    # Tune the proposal covariance for each transition to achieve
                    # reasonable acceptance ratios.
                    if(ttt %% 30 == 0){
                        if(ttt %% 480 == 0){
                            accept[j] = 0
                            
                        } else if( accept[j] / (ttt %% 480) < .4 ){ 
                            pscale[j] = (.75^2)*pscale[j]
                            
                        } else if( accept[j] / (ttt %% 480) > .5 ){ 
                            pscale[j] = (1.25^2)*pscale[j]
                        }
                    }
                }
                # --------------------------------------------------------------
        
            } else {
                # MH update: R -------------------------------------------------
        
                # Prior for R 
                nu_R = 1000
                psi_R = diag(c(9, 81, 81, 9))
                psi_R = (nu_R - 4 - 1) * psi_R
        
                # q(R* | R(t)) 
                curr_R = matrix(par[ind_j], nrow = 4)
                psi_nu_star_t = proposal_R_cpp_new(nu_R, psi_R, curr_R, Y, Dn, Xn,
                                                   A, par, par_index, as.numeric(EIDs),
                                                   B, Dn_omega, W)
                q_s_star_t = psi_nu_star_t[[1]] / pscale[j]
                q_nu_star_t = floor(psi_nu_star_t[[2]] / pscale[j])
        
                # Proposal R 
                proposal[ind_j] = c(rinvwishart(nu = q_nu_star_t, S = q_s_star_t))
        
                # q(R(t) | R*) 
                prop_R = matrix(proposal[ind_j], nrow = 4)
                psi_nu_t_star = proposal_R_cpp_new(nu_R, psi_R, prop_R, Y, Dn, Xn,
                                                   A, par, par_index, as.numeric(EIDs),
                                                   B, Dn_omega, W)
                q_s_t_star = psi_nu_t_star[[1]] / pscale[j]
                q_nu_t_star = floor(psi_nu_t_star[[2]] / pscale[j])
        
                # Evaluate log-density and log-likelihood
                log_prop      = dinvwishart(Sigma = prop_R,
                                            nu = q_nu_star_t,
                                            S = q_s_star_t, log = T)
                log_prop_prev = dinvwishart(Sigma = curr_R,
                                            nu = q_nu_t_star,
                                            S = q_s_t_star, log = T)
        
                log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index,
                                           A, B, Y, z, Dn, Xn, Dn_omega, W, n_cores)
        
                # Catch any numerical issues from C++
                if(!is.finite(log_target + log_prop_prev - log_target_prev - log_prop) | 
                    is.nan(log_target + log_prop_prev - log_target_prev - log_prop)) {
                    print("Current R")
                    print(curr_R)
                    print("Prop R")
                    print(prop_R)
                    print(paste0("log_target: ", log_target))
                    print(paste0("log_prop_prev: ", log_prop_prev))
                    print(paste0("log_target_prev: ", log_target_prev))
                    print(paste0("log_prop: ", log_prop))
                    print(log_target + log_prop_prev - log_target_prev - log_prop)
                }
        
                # Compute Metropolis ratio 
                if( log_target + log_prop_prev - log_target_prev - log_prop > log(runif(1,0,1)) ){
                    log_target_prev = log_target
                    par[ind_j] = proposal[ind_j]
                    accept[j] = accept[j] +1
                }
        
                chain[chain_ind,ind_j] = par[ind_j]
            }
        }

        # Restart the acceptance ratio at burnin
        if(ttt == burnin) accept = rep( 0, n_group)
        
        B_chain[ chain_ind, ] = do.call( 'c', B)
        hc_chain[chain_ind, ] = Y[,'hemo']
        hr_chain[chain_ind, ] = Y[,'hr']
        bp_chain[chain_ind, ] = Y[,'map']
        la_chain[chain_ind, ] = Y[,'lactate']
        # ----------------------------------------------------------------------

        if(ttt%%1==0)  cat('--->',ttt,'\n')
        if(ttt > burnin & ttt%%chain_length_MASTER==0) {
            mcmc_end_t = Sys.time() - mcmc_start_t; print(mcmc_end_t)
            index_keep = seq(1, chain_length_MASTER, by = 5)
            mcmc_out_temp = list(chain    = chain[index_keep,], 
                                 B_chain  = B_chain[index_keep,], 
                                 hc_chain = hc_chain[index_keep,], 
                                 hr_chain = hr_chain[index_keep,],
                                 bp_chain = bp_chain[index_keep,], 
                                 la_chain = la_chain[index_keep,],
                                 # A_chain  = A_chain,
                                 otype=otype, accept=accept/length(burnin:ttt), 
                                 pscale=pscale, pcov = pcov, par_index=par_index)
            if(simulation) {
                save(mcmc_out_temp, file = paste0('Model_out/mcmc_out_',ind,'_', 
                                                  trialNum, 'it', ttt/chain_length_MASTER + (max_ind - 5), 
                                                  '_samp', sampling_num, '_sim.rda'))
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



# Create function for initializing missing Y
# if(max_ind == 5) {
#     file_name = paste0("Data_updates/Y_init", trialNum, "_df", df_num, ".rda")
#     if(file.exists(file_name)) {
#         load(file_name)
#     } else {
#         # Setting initial values for Y
#         print("Initializing the missing Y's for imputation")
#         for(i in EIDs) {
#             print(which(EIDs == i))
#             heading_names = c('hemo', 'hr', 'map', 'lactate')
#             noise_var = c(diag(matrix(par[par_index$vec_R], nrow = 4)))
#             sub_dat = Y[Y[,"EID"] == i, ]
    
#             for(k in 1:length(heading_names)) {
#                 if(sum(is.na(sub_dat[,heading_names[k]])) == nrow(sub_dat)) {
#                     sub_dat[,heading_names[k]] = mean(Y[,heading_names[k]], na.rm =T) + 
#                         rnorm(n = nrow(sub_dat), mean = 0, sd = sqrt(noise_var[k]))
#                 } else {
#                     if(sum(!is.na(sub_dat[,heading_names[k]])) == 1) {
#                         non_na_val = sub_dat[!is.na(sub_dat[,heading_names[k]]), heading_names[k]]
#                         sub_dat[is.na(sub_dat[,heading_names[k]]), heading_names[k]] = rnorm(n = sum(is.na(sub_dat[,heading_names[k]])), 
#                                                                                                 mean = non_na_val,
#                                                                                                 sd = sqrt(noise_var[k]))
#                         # sub_dat[,heading_names[k]] = sub_dat[!is.na(sub_dat[,heading_names[k]]), heading_names[k]]
#                     } else {
#                         obs_indices = which(!is.na(sub_dat[,heading_names[k]]))
#                         miss_indices = which(is.na(sub_dat[,heading_names[k]]))
#                         for(j in miss_indices) {
#                             if(j < obs_indices[1]) {
#                                 sub_dat[j,heading_names[k]] = sub_dat[obs_indices[1], heading_names[k]] + 
#                                     rnorm(n = 1, mean = 0, sd = sqrt(noise_var[k]))
#                             } else if(j > tail(obs_indices,1)) {
#                                 sub_dat[j,heading_names[k]] = sub_dat[tail(obs_indices,1), heading_names[k]] + 
#                                     rnorm(n = 1, mean = 0, sd = sqrt(noise_var[k]))
#                             } else {
#                                 end_pts = c(max(obs_indices[obs_indices < j]),
#                                             min(obs_indices[obs_indices > j]))
#                                 slope = (sub_dat[end_pts[2], heading_names[k]] - sub_dat[end_pts[1], heading_names[k]]) / diff(end_pts)
#                                 sub_dat[j,heading_names[k]] = slope * (j - end_pts[1]) + sub_dat[end_pts[1], heading_names[k]] + 
#                                                                         rnorm(n = 1, mean = 0, sd = sqrt(noise_var[k]))
#                             }
#                         }
#                     }
#                 }
#                 Y[Y[,"EID"] == i, heading_names[k]] = sub_dat[,heading_names[k]]
#             }
#         }
        
#         # Ensure we don't have any negative values
#         Y[Y[,'hemo']    < 0,'hemo']    = 0.001
#         Y[Y[,'hr']      < 0,'hr']      = 0.001
#         Y[Y[,'map']     < 0,'map']     = 0.001
#         Y[Y[,'lactate'] < 0,'lactate'] = 0.001
        
#         save(Y, file = file_name)
#         print("Done initializing")
#     }
# } else {
#     print("continue from last iteration")
#     # ----------------------------------------------------------------------
#     prev_file = paste0('Model_out/mcmc_out_interm_', ind, '_', trialNum, 'it', max_ind-5, '_df', df_num, '.rda')
#     load(prev_file)
#     # ----------------------------------------------------------------------
    
#     pcov = mcmc_out_temp$pcov
#     pscale = mcmc_out_temp$pscale
    
#     # Setting initial values for Y
#     Y[, 'hemo']    = c(mcmc_out_temp$hc_chain[nrow(mcmc_out_temp$hc_chain), ])
#     Y[, 'hr']      = c(mcmc_out_temp$hr_chain[nrow(mcmc_out_temp$hr_chain), ])
#     Y[, 'map']     = c(mcmc_out_temp$bp_chain[nrow(mcmc_out_temp$bp_chain), ])
#     Y[, 'lactate'] = c(mcmc_out_temp$la_chain[nrow(mcmc_out_temp$la_chain), ])
    
#     # Ensure we don't have any negative values
#     Y[Y[,'hemo']    < 0,'hemo']    = 0.001
#     Y[Y[,'hr']      < 0,'hr']      = 0.001
#     Y[Y[,'map']     < 0,'map']     = 0.001
#     Y[Y[,'lactate'] < 0,'lactate'] = 0.001
    
#     rm(mcmc_out_temp)
# }


# Initialize B based on max likelihood at each point given initial values --
# if(max_ind == 5) {
#     Xn_initial = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x, 10)
#     Xn = Xn_initial[[2]]
#     
#     B = list()
#     B = initialize_b_i(as.numeric(EIDs), par, par_index, A, Y, z, Xn, Dn_omega, W, 10)
#     
#     # Reset the clinic_rule subjects
#     if(!simulation) {
#         c_r = unique(Y[Y[,"clinic_rule"] < 0, "EID"])  
#         for(j in c_r) {
#             c_r_ind = which(EIDs == j)
#             B[[c_r_ind]] = matrix(1, nrow = sum(Y[,"EID"] == j), ncol = 1)
#         }
#     }
# }
# --------------------------------------------------------------------------