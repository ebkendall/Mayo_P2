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
                        par, par_index, Y, B, A, y_first){
    
    EIDs = as.numeric(unique(Y[,'EID']))
    
    # Number of cores over which to parallelize --------------------------------
    n_cores = strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
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

    # Metropolis Parameter Index for MH within Gibbs updates -------------------
    mpi = list(c(par_index$init),
               c(par_index$zeta),
               c(par_index$A),
               c(par_index$R))

    n_group = length(mpi)
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))
    pscale = rep( 0.001, n_group)
  
    # Initialize data storage --------------------------------------------------
    reset_step = 10000
    chain_length_MASTER = 1000

    chain = matrix(NA, reset_step, length(par)) 
    B_chain = matrix(NA, chain_length_MASTER, nrow(Y)) 

    accept = rep( 0, n_group)
    
    # Initialization continued -------------------------------------------------
    B_Dn = mle_state_seq(EIDs, par, par_index, A, Y, y_first, n_cores)
    B = B_Dn[[1]]
    Dn = B_Dn[[2]]
    # Dn = initialize_Dn(EIDs, B)

    # Start Metropolis-within-Gibbs Algorithm ----------------------------------
    chain[1,] = par
    B_chain[1, ] = do.call( 'c', B)
    
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
        gamma_1 = update_gamma_i(EIDs, par, par_index, A, B, Dn, Y, n_cores, y_first)

        # State-space update (B) -----------------------------------------------
        sps = sample(x = 10:50, size = 1, replace = T) # sps >= 2
        B_Dn = state_sampler(EIDs, par, par_index, A, B, Y, sps, gamma_1, n_cores)
        B = B_Dn[[1]]
        Dn = B_Dn[[2]]
        
        # Gibbs: alpha_i -------------------------------------------------------
        A = update_alpha_i(EIDs, par, par_index, B, Dn, Y, n_cores, y_first, gamma_1)
        
        # Gibbs: alpha~, omega~, beta, & Upsilon -------------------------------
        par = update_alpha_tilde(EIDs, par, par_index, A, Y)
        par = update_beta_upsilon(EIDs, par, par_index, A, n_cores)

        # Store current parameter updates --------------------------------------
        chain[chain_ttt,] = par
        
        # Evaluate log-likelihood before MH step -------------------------------
        log_target_prev = log_post(EIDs, par, par_index, A, B, Dn, Y, n_cores, 
                                   y_first, gamma_1)

        if(!is.finite(log_target_prev)){
            print(paste0("Infinite log-posterior: ", log_target_prev)); stop();
        }

        # Metropolis-Hastings updates ------------------------------------------
        for(j in 1:n_group) {

            ind_j = mpi[[j]]
            proposal = par

            # MH update: A, zeta, init -----------------------------------------
            if(sum(ind_j %in% par_index$R) == 0) {
                proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                           sigma=pscale[[j]]*pcov[[j]])

                # Evaluate proposed log-likelihood -----------------------------
                log_target = log_post(EIDs, proposal, par_index, A, B, Dn, Y, 
                                      n_cores, y_first, gamma_1)

                if(ttt < burnin){
                    while(!is.finite(log_target)){
                        print('bad proposal')
                        proposal = par
                        proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                                   sigma=pcov[[j]]*pscale[j])

                        log_target = log_post(EIDs, proposal, par_index, A, B, Dn, Y, 
                                              n_cores, y_first, gamma_1)
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
                # --------------------------------------------------------------

            } else {
                # MH update: R -------------------------------------------------
                curr_R = matrix(par[ind_j], nrow = 4)
                
                # Prior for R
                nu_R = 8
                psi_R = diag(c(9, 9, 9, 9))
                psi_R = (nu_R - 4 - 1) * psi_R
                
                # Proposal
                psi_nu_q_star = proposal_R(nu_R, psi_R, curr_R, EIDs, par, 
                                           par_index, A, B, Dn, Y, n_cores, 
                                           y_first, gamma_1)
                psi_q_star = psi_nu_q_star[[1]]
                nu_q_star = psi_nu_q_star[[2]]
                
                proposal[ind_j] = c(rinvwishart(nu = nu_q_star, S = psi_q_star))
                prop_R = matrix(proposal[ind_j], nrow = 4)
                
                # Current
                psi_nu_q_curr = proposal_R(nu_R, psi_R, prop_R, EIDs, par, 
                                           par_index, A, B, Dn, Y, n_cores, 
                                           y_first, gamma_1)
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
                log_target = log_post(EIDs, proposal, par_index, A, B, Dn, Y, 
                                      n_cores, y_first, gamma_1)

                if(!is.finite(log_target + log_q_curr - log_target_prev - log_q_prop) |
                   is.nan(log_target + log_q_curr - log_target_prev - log_q_prop)) {
                    print("Current R")
                    print(curr_R)
                    print("Prop R")
                    print(prop_R)
                    print(paste0("log_target: ", log_target))
                    print(paste0("log_q_curr: ", log_q_curr))
                    print(paste0("log_target_prev: ", log_target_prev))
                    print(paste0("log_q_prop: ", log_q_prop))
                    print(log_target + log_q_curr - log_target_prev - log_q_prop)
                }

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
            # index_keep_2 = seq(2, chain_length_MASTER, by = 2)
            
            if(simulation) {
                mcmc_out = list(chain    = chain[index_keep,], 
                                B_chain  = B_chain, 
                                hc_chain = Y[,'hemo'],
                                hr_chain = Y[,'hr'], 
                                bp_chain = Y[,'map'], 
                                la_chain = Y[,'lactate'], 
                                accept=accept/length((burnin+1):ttt), 
                                pscale=pscale, pcov = pcov, par_index=par_index)

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
                                A_chain  = A_chain,
                                W_chain  = W_chain,
                                otype=otype, accept=accept/length((burnin+1):ttt), 
                                pscale=pscale, pcov = pcov, par_index=par_index)

                save(mcmc_out, 
                     file = paste0('Model_out/mcmc_out_', trialNum, '_', seed_num, 
                                   'it', ttt/reset_step + (max_ind - 5), '.rda'))
            }
            
            # Reset the chains
            chain = matrix(NA, reset_step, length(par)) 
            B_chain = matrix(NA, chain_length_MASTER, nrow(Y)) 
        }
    }
    # --------------------------------------------------------------------------

    return(0)
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------