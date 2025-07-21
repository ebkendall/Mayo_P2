library(mvtnorm, quietly=T)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(expm, quietly = T)
sourceCpp("mcmc_fnc.cpp")

mcmc_routine = function(par, par_index, B, y, ids, steps, burnin, ind, dgm){
    
    EIDs = unique(ids)
    
    # Number of cores over which to parallelize --------------------------------
    n_cores = 6
    print(paste0("Number of cores: ", n_cores))
    
    # Transition information ---------------------------------------------------
    adjacency_mat = matrix(c(1,1,0,
                             0,1,1,
                             1,1,1), ncol = 3, byrow = T)
    initialize_cpp(adjacency_mat)
    
    # Metropolis Parameter Index for MH within Gibbs updates -------------------
    if(dgm) {
        mpi = list(c(par_index$alpha[c(1,4,7,10)]), 
                   c(par_index$alpha[c(2,5,8,11)]),
                   c(par_index$alpha[c(3,6,9,12)]),
                   c(par_index$zeta),
                   c(par_index$diag_R),
                   c(par_index$init))
    } else {
        mpi = list(c(par_index$alpha[c(1,3,5,7)]), 
                   c(par_index$alpha[c(2,4,6,8)]),
                   c(par_index$zeta),
                   c(par_index$diag_R),
                   c(par_index$init),
                   c(par_index$diag_G))
    }
    
    n_group = length(mpi)
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))
    pscale = rep(0.001, n_group)
    
    # Initialize data storage --------------------------------------------------
    reset_step = 10000
    chain_length_MASTER = 1000
    
    chain = matrix(NA, reset_step, length(par)) 
    B_chain = matrix(NA, chain_length_MASTER, nrow(y)) 
    
    accept = rep( 0, n_group)
    
    # Initialize states to MLE state sequences ---------------------------------
    B = mle_state_seq(as.numeric(EIDs), par, par_index, y, ids, dgm)
    
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
        
        # Sample gamma_i's -----------------------------------------------------
        gamma_i = matrix(0, nrow = 1, ncol = 4)
        if(!dgm) {
            gamma_i = gamma_i_sample(as.numeric(EIDs), par, par_index, B, y, ids, n_cores)
        }

        # Efficient state-sampler ----------------------------------------------
        sps = sample(x = 2:50, size = 1, replace = T) # sps >= 2
        B_Dn = state_sampler(as.numeric(EIDs), par, par_index, B, y, ids, n_cores, sps, dgm, gamma_i)
        B = B_Dn

        # Evaluate log-likelihood before MH step -------------------------------
        log_target_prev = log_post_cpp(as.numeric(EIDs), par, par_index, B,
                                       y, ids, n_cores, dgm, gamma_i)

        if(!is.finite(log_target_prev)){
            print(paste0("Infinite log-posterior: ", log_target_prev)); stop();
        }

        # Metropolis-Hastings updates ------------------------------------------
        for(j in 1:n_group) {

            ind_j = mpi[[j]]
            proposal = par

            if(length(ind_j) > 1) {
                proposal[ind_j] = rmvnorm(n=1, mean=par[ind_j],
                                       sigma=pscale[[j]]*pcov[[j]])
            } else {
                proposal[ind_j] = rnorm(n=1, mean=par[ind_j],
                                        sd=sqrt(pscale[[j]]*pcov[[j]]))
            }

            # Evaluate proposed log-likelihood -----------------------------
            log_target = log_post_cpp(as.numeric(EIDs), proposal, par_index,
                                      B, y, ids, n_cores, dgm, gamma_i)

            if(ttt < burnin){
                while(!is.finite(log_target)){
                    print('bad proposal')
                    proposal = par

                    if(length(ind_j) > 1) {
                        proposal[ind_j] = rmvnorm(n=1, mean=par[ind_j],
                                            sigma=pscale[[j]]*pcov[[j]])
                    } else {
                        proposal[ind_j] = rnorm(n=1, mean=par[ind_j],
                                                sd=sqrt(pscale[[j]]*pcov[[j]]))
                    }

                    log_target = log_post_cpp(as.numeric(EIDs), proposal, par_index, 
                                              B, y, ids, n_cores, dgm, gamma_i)
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

            chain[chain_ttt, ind_j] = par[ind_j]

            # Proposal tuning scheme ---------------------------------------
            # During the burnin period, update the proposal covariance in each step
            # to capture the relationships within the parameters vectors for each
            # transition.  This helps with mixing.
            if(ttt < burnin){
                if(ttt == 100)  pscale[j] = 1

                if(100 <= ttt & ttt <= 2000){
                    temp_chain = chain[1:ttt,ind_j]
                    if(length(ind_j) > 1) {
                        pcov[[j]] = cov(temp_chain[!duplicated(temp_chain),, drop=F])
                    } else {
                        pcov[[j]] = var(temp_chain[!duplicated(temp_chain)])
                    }

                } else if(2000 < ttt){
                    temp_chain = chain[(ttt-2000):ttt,ind_j]
                    if(length(ind_j) > 1) {
                        pcov[[j]] = cov(temp_chain[!duplicated(temp_chain),, drop=F])
                    } else {
                        pcov[[j]] = var(temp_chain[!duplicated(temp_chain)])
                    }
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
        }

        chain[chain_ttt, ] = par
        
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
        
        if(ttt > burnin & ttt%%reset_step == 0) {
            mcmc_end_t = Sys.time() - mcmc_start_t; print(mcmc_end_t)
            
            index_keep = seq(10, reset_step, by = 10)
            
            mcmc_out = list( chain    = chain[index_keep,], 
                             B_chain  = B_chain, 
                             accept=accept/length(burnin:ttt), 
                             pscale=pscale, pcov = pcov, par_index=par_index)
            
            save(mcmc_out, file = paste0('Model_out/mcmc_out_', ind, '_it_', ttt/reset_step,
                                         '_', as.numeric(dgm), '.rda'))
            
            # Reset the chains
            chain = matrix(NA, reset_step, length(par)) 
            B_chain = matrix(NA, chain_length_MASTER, nrow(y)) 
        }
    }
    # --------------------------------------------------------------------------
    
    return(list( chain=chain, B_chain=B_chain, accept=accept/(steps-burnin), 
                 pscale=pscale, pcov = pcov, par_index=par_index))
}
