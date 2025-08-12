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
                        par, par_index, Y, x, z, B, A, Dn_omega, bleed_indicator){
    
    EIDs = as.numeric(unique(Y[,'EID']))
    
    # Number of cores over which to parallelize --------------------------------
    n_cores = 6
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
    
    # Metropolis Parameter Index for MH within Gibbs updates -------------------
    mpi = list(c(par_index$init),
               c(par_index$zeta[1:2]),   c(par_index$zeta[3:4]),
               c(par_index$zeta[5:6]),   c(par_index$zeta[7:8]),
               c(par_index$zeta[9:10]),  c(par_index$zeta[11:12]),
               c(par_index$zeta[13:14]), c(par_index$zeta[15:16]),
               c(par_index$zeta[17:18]), c(par_index$zeta[19:20]),
               c(par_index$zeta[21:22]), c(par_index$zeta[23:24]),
               c(par_index$A[1:4]), c(par_index$A[5:8]), c(par_index$A[9:12]), 
               c(par_index$A[13:16]), c(par_index$A[17:20]),
               c(par_index$R))
    
    n_group = length(mpi)
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))
    pscale = rep( 0.001, n_group)
    accept = rep( 0, n_group)
    
    # Initialize design matrix -------------------------------------------------
    Xn = initialize_Xn(EIDs, Y, x)
    
    # Initialize Y and states --------------------------------------------------
    impute_step = TRUE
    if(sum(c(otype) == 0) > 0) {
        
        # There exists missing Y values
        print("Missing Y values - initializing")
        
        if(max_ind > 5) {
            
            chosen_seed = 10
            
            load(paste0('Model_out/mcmc_out_', trialNum, '_', chosen_seed, 'it', 
                        max_ind - 5, '.rda'))
            
            # initialize Y
            Y[,'hemo'] = mcmc_out$hc_chain[nrow(mcmc_out$hc_chain), ]
            Y[,'hr'] = mcmc_out$hr_chain[nrow(mcmc_out$hr_chain), ]
            Y[,'map'] = mcmc_out$bp_chain[nrow(mcmc_out$bp_chain), ]
            Y[,'lactate'] = mcmc_out$la_chain[nrow(mcmc_out$la_chain), ]
            
            # initialize proposal structure
            pcov = mcmc_out$pcov
            pscale = mcmc_out$pscale
            
            # initialize Dn_alpha
            Dn_alpha = initialize_Dn(EIDs, B)
            
            rm(mcmc_out)
            
        } else {
            
            vital_means = colMeans(Y[,c('hemo', 'hr', 'map', 'lactate')], na.rm = T)
            Y_init = initialize_Y(EIDs, par, par_index, A, Y, z, Dn_omega, 
                                  Xn, otype, n_cores, vital_means)
            
            Y[, c('hemo', 'hr', 'map', 'lactate')] = Y_init
            
            B_Dn = mle_state_seq(EIDs, par, par_index, A, Y, z, Dn_omega, Xn, n_cores)
            B = B_Dn[[1]]
            Dn_alpha = B_Dn[[2]]
        }
        
        print("Done")
        
    } else {
        
        # There are NO missing Y values
        print("No missingness")
        impute_step = FALSE
        
        B_Dn = mle_state_seq(EIDs, par, par_index, A, Y, z, Dn_omega, Xn, n_cores)
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
        gamma_1 = update_gamma_i(EIDs, par, par_index, A, B, Y, Dn_alpha, Dn_omega, Xn, n_cores)
        
        # Imputing the missing Y values ----------------------------------------
        if(impute_step) {
            Y = impute_Y(EIDs, par, par_index, A, B, Y, Dn_alpha, Dn_omega, Xn, 
                         gamma_1, otype, n_cores)
            colnames(Y) = c('EID','hemo', 'hr', 'map', 'lactate',
                            'RBC_rule', 'clinic_rule')    
        }
        
        # # State-space update (B) -----------------------------------------------
        # sps = sample(x = 2:50, size = 1, replace = T) # sps >= 2
        # B_Dn = state_sampler(EIDs, par, par_index, B, A, Y, z, Dn_alpha, 
        #                      Dn_omega, Xn, bleed_indicator, gamma_1, sps, n_cores)
        # B = B_Dn[[1]]
        # Dn_alpha = B_Dn[[2]]
        
        # Gibbs: alpha_i -------------------------------------------------------
        A = update_alpha_i(EIDs, par, par_index, B, Y, Dn_alpha, Dn_omega, Xn, gamma_1, n_cores)
        
        # Gibbs: alpha~, omega~, beta, Upsilon, and G --------------------------
        par = update_alpha_tilde(EIDs, par, par_index, A, Y)
        par = update_omega(EIDs, par, par_index, A, B, Y, Dn_alpha, Dn_omega, 
                           Xn, gamma_1, n_cores)
        par = update_beta_upsilon(EIDs, par, par_index, A, B, Y, Dn_alpha, Dn_omega,
                                  Xn, gamma_1, n_cores)
        par = update_G(EIDs, par, par_index, Y, Dn_omega, Xn, gamma_1, n_cores)
        
        # Store current parameter updates --------------------------------------
        chain[chain_ttt,] = par
        
        # Evaluate log-likelihood before MH step -------------------------------
        log_target_prev = log_post(EIDs, par, par_index, B, A, Y, z, Dn_alpha, 
                                   Dn_omega, Xn, gamma_1, n_cores)
        
        if(!is.finite(log_target_prev)){
            print(paste0("Infinite log-posterior: ", log_target_prev)); stop();
        }
        
        # Metropolis-Hastings updates ------------------------------------------
        for(j in 1:n_group) {
            
            ind_j = mpi[[j]]
            proposal = par
            
            # MH update: init, zeta, A, eta ------------------------------------
            if(sum(ind_j %in% par_index$R) == 0) {
                
                proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                           sigma=pscale[[j]]*pcov[[j]])
                
                # Evaluate proposed log-likelihood -----------------------------
                log_target = log_post(EIDs, proposal, par_index, B, A, Y, z, Dn_alpha, 
                                      Dn_omega, Xn, gamma_1, n_cores)
                
                if(ttt < burnin){
                    while(!is.finite(log_target)){
                        print('bad proposal')
                        proposal = par
                        proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                                   sigma=pcov[[j]]*pscale[j])
                        
                        log_target = log_post(EIDs, proposal, par_index, B, A, Y, z, Dn_alpha, 
                                              Dn_omega, Xn, gamma_1, n_cores)
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
                
                chain[chain_ttt,ind_j] = par[ind_j]
                
                # Proposal tuning scheme ---------------------------------------
                # During the burnin period, update the proposal covariance
                # to capture the relationships within the parameter vectors.
                # This helps with mixing.
                if(ttt < burnin){
                    if(chain_ttt == 100)  pscale[j] = 1
                    
                    if(100 <= chain_ttt & chain_ttt <= 2000){
                        temp_chain = chain[1:chain_ttt,ind_j]
                        pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                        
                    } else if(2000 < chain_ttt){
                        temp_chain = chain[(chain_ttt-2000):chain_ttt,ind_j]
                        pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                    }
                    if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )
                    
                    # Tune the proposal covariance for each transition to achieve
                    # reasonable acceptance ratios.
                    if(chain_ttt %% 30 == 0){
                        if(chain_ttt %% 480 == 0){
                            accept[j] = 0
                            
                        } else if( accept[j] / (chain_ttt %% 480) < .4 ){
                            pscale[j] = (.75^2)*pscale[j]
                            
                        } else if( accept[j] / (chain_ttt %% 480) > .5 ){
                            pscale[j] = (1.25^2)*pscale[j]
                        }
                    }
                }
            } else {
                
                # MH update: R -------------------------------------------------
                curr_R = matrix(par[ind_j], nrow = 4)
                
                # Prior for R
                nu_R = 8
                psi_R = diag(c(4, 16, 16, 4))
                psi_R = (nu_R - 4 - 1) * psi_R
                
                # Proposal
                psi_nu_q_star = proposal_R(nu_R, psi_R, curr_R, EIDs, par, 
                                           par_index, A, B, Y, Dn_alpha, 
                                           Dn_omega, Xn, gamma_1, n_cores)
                psi_q_star = psi_nu_q_star[[1]]
                nu_q_star = psi_nu_q_star[[2]]
                
                proposal[ind_j] = c(rinvwishart(nu = nu_q_star, S = psi_q_star))
                prop_R = matrix(proposal[ind_j], nrow = 4)
                
                # Current
                psi_nu_q_curr = proposal_R(nu_R, psi_R, prop_R, EIDs, par, 
                                           par_index, A, B, Y, Dn_alpha, 
                                           Dn_omega, Xn, gamma_1, n_cores)
                psi_q_curr = psi_nu_q_curr[[1]]
                nu_q_curr = psi_nu_q_curr[[2]]
                
                # q(R* | R)
                log_q_prop = dinvwishart(Sigma = prop_R, 
                                         nu = nu_q_star, 
                                         S = psi_q_star, log = T)
                
                # q(R | R*)
                log_q_curr = dinvwishart(Sigma = curr_R, 
                                         nu = nu_q_curr, 
                                         S = psi_q_curr, log = T)
                
                # Log-posterior at proposal
                log_target = log_post(EIDs, proposal, par_index, B, A, Y, z, Dn_alpha, 
                                      Dn_omega, Xn, gamma_1, n_cores)
                
                # Compute Metropolis ratio
                if( log_target + log_q_curr - log_target_prev - log_q_prop > log(runif(1,0,1)) ){
                    log_target_prev = log_target
                    par[ind_j] = proposal[ind_j]
                    accept[j] = accept[j] +1
                }
                
                chain[chain_ttt,ind_j] = par[ind_j]
            }
        }
        
        # Restart the acceptance ratio at burnin
        if(ttt == burnin) accept = rep( 0, n_group)
        
        B_chain[ chain_ind, ] = do.call( 'c', B)
        if(impute_step) {
            hc_chain[chain_ind, ] = Y[,'hemo']
            hr_chain[chain_ind, ] = Y[,'hr']
            bp_chain[chain_ind, ] = Y[,'map']
            la_chain[chain_ind, ] = Y[,'lactate']
        }
        # ----------------------------------------------------------------------
        
        cat('--> ', ttt, ', c_ttt = ', chain_ttt, ', c_ind = ', chain_ind, '\n')
        
        if(ttt%%100==0) {
            print(accept) 
            print(pscale)
        }
        
        ttt_end_t = Sys.time() - ttt_start_t; print(ttt_end_t)
        
        if(ttt > burnin & ttt%%reset_step==0) {
            mcmc_end_t = Sys.time() - mcmc_start_t; print(mcmc_end_t)
            
            index_keep = seq(10, reset_step, by = 10)
            index_keep_2 = seq(2, chain_length_MASTER, by = 2)
            
            if(simulation) {
                if(impute_step) {
                    mcmc_out = list(chain    = chain[index_keep,], 
                                    B_chain  = B_chain[index_keep_2,], 
                                    hc_chain = hc_chain[index_keep_2,], 
                                    hr_chain = hr_chain[index_keep_2,], 
                                    bp_chain = bp_chain[index_keep_2,], 
                                    la_chain = la_chain[index_keep_2,], 
                                    alpha_i = A,
                                    otype=otype, accept=accept/length((burnin+1):ttt), 
                                    pscale=pscale, pcov = pcov, par_index=par_index)    
                } else {
                    mcmc_out = list(chain    = chain[index_keep,], 
                                    B_chain  = B_chain[index_keep_2,], 
                                    hc_chain = matrix(Y[,'hemo'], nrow = 1),
                                    hr_chain = matrix(Y[,'hr'], nrow = 1), 
                                    bp_chain = matrix(Y[,'map'], nrow = 1), 
                                    la_chain = matrix(Y[,'lactate'], nrow = 1), 
                                    alpha_i = A,
                                    otype=otype, accept=accept/length((burnin+1):ttt), 
                                    pscale=pscale, pcov = pcov, par_index=par_index)
                }
                
                save(mcmc_out, 
                     file = paste0('Model_out/mcmc_out_', trialNum, '_', seed_num,
                                   'it', ttt/reset_step + (max_ind - 5), '_sim.rda'))
            } else {
                mcmc_out = list(chain    = chain[index_keep,], 
                                B_chain  = B_chain[index_keep_2,], 
                                hc_chain = hc_chain[index_keep_2,], 
                                hr_chain = hr_chain[index_keep_2,],
                                bp_chain = bp_chain[index_keep_2,], 
                                la_chain = la_chain[index_keep_2,],
                                alpha_i = A,
                                otype=otype, accept=accept/length((burnin+1):ttt), 
                                pscale=pscale, pcov = pcov, par_index=par_index)
                
                save(mcmc_out, 
                     file = paste0('Model_out/mcmc_out_', trialNum,'_', seed_num, 
                                   'it', ttt/reset_step + (max_ind - 5), '.rda'))
            }
            
            # Reset the chains
            chain = matrix(NA, reset_step, length(par)) 
            B_chain = matrix(NA, chain_length_MASTER, nrow(Y)) 
            if(impute_step) {
                hc_chain = hr_chain = matrix(NA, chain_length_MASTER, nrow(Y))
                bp_chain = la_chain = matrix(NA, chain_length_MASTER, nrow(Y))    
            }
        }
    }
    # ---------------------------------------------------------------------------
    
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
