library(mvtnorm, quietly=T)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(expm, quietly = T)
sourceCpp("mcmc_fnc.cpp")

mcmc_routine = function(par, par_index, B, y, ids, steps, burnin, ind, 
                        sampling_num, states_per_step, steps_per_it){
    
    EIDs = unique(ids)
    
    # Number of cores over which to parallelize --------------------------------
    n_cores = 8#strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
    print(paste0("Number of cores: ", n_cores))
    
    # Transition information ---------------------------------------------------
    adjacency_mat = matrix(c(1,1,
                             1,1), ncol = 2, byrow = T)
    initialize_cpp(adjacency_mat, states_per_step)
    
    # Metropolis Parameter Index for MH within Gibbs updates -------------------
    mpi = list(c(par_index$mu), 
               c(par_index$t_p))
    
    n_group = length(mpi)
    pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(mpi[[j]]))
    pscale = rep( 0.0001, n_group)
    
    # Initialize data storage --------------------------------------------------
    chain_length_MASTER = 10000
    chain = matrix( 0, chain_length_MASTER, length(par)) 
    B_chain = matrix( 0, chain_length_MASTER, length(y)) 
    accept = rep( 0, n_group)
    
    # Initialize states to MLE state sequences ---------------------------------
    B = mle_state_seq(as.numeric(EIDs), par, par_index, y, ids, n_cores)
    
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
        
        for(s in 1:steps_per_it) {
            # Random sample update ---------------------------------------------
            if(sampling_num == 1) {
                B_Dn = mh_up(as.numeric(EIDs), par, par_index, B, y, ids,
                             n_cores, states_per_step)
                B = B_Dn
            }

            # Almost-Gibbs update ----------------------------------------------
            if(sampling_num == 2) {
                B_Dn = almost_gibbs_up(as.numeric(EIDs), par, par_index, B, y,
                                       ids, n_cores, states_per_step)
                B = B_Dn
            }

            # Gibbs update -----------------------------------------------------
            if(sampling_num == 3) {
                B_Dn = gibbs_up(as.numeric(EIDs), par, par_index, B, y, ids,
                                n_cores, states_per_step)
                B = B_Dn
            }
            
            # Full seq MH update -----------------------------------------------
            if(sampling_num == 4) {
                B_Dn = mh_up_all(as.numeric(EIDs), par, par_index, B, y, ids,
                                 n_cores)
                B = B_Dn
            }
        }

        # Evaluate log-likelihood before MH step -------------------------------
        log_target_prev = log_post_cpp(as.numeric(EIDs), par, par_index, B,
                                       y, ids, n_cores)

        if(!is.finite(log_target_prev)){
            print("Infinite log-posterior")
            print(paste0("value: ", log_target_prev))
            break
        }

        # Metropolis-Hastings updates ------------------------------------------
        for(j in 1:n_group) {

            ind_j = mpi[[j]]
            proposal = par

            proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                       sigma=pscale[[j]]*pcov[[j]])

            # Evaluate proposed log-likelihood -----------------------------
            log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, B, y, ids, n_cores)

            if(ttt < burnin){
                while(!is.finite(log_target)){
                    print('bad proposal')
                    proposal = par
                    proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                               sigma=pcov[[j]]*pscale[j])

                    log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, B, y, ids, n_cores)
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

                    } else if( accept[j] / (ttt %% 480) < .45 ){
                        pscale[j] = (.5^2)*pscale[j]

                    } else if( accept[j] / (ttt %% 480) > .55 ){
                        pscale[j] = (1.5^2)*pscale[j]
                    }
                }
            }
            # --------------------------------------------------------------
        }
        
        # Restart the acceptance ratio at burnin
        if(ttt == burnin) accept = rep( 0, n_group)
        
        B_chain[ chain_ind, ] = do.call( 'c', B)
        # ----------------------------------------------------------------------
        
        if(ttt%%1==0)  print(paste0('--->',ttt))
        if(ttt%%100==0) print(accept)
        
        if(ttt > burnin & ttt%%chain_length_MASTER==0) {
            mcmc_end_t = Sys.time() - mcmc_start_t; print(mcmc_end_t)
            index_keep = seq(1, chain_length_MASTER, by = 5)
            mcmc_out = list( chain    = chain[index_keep,], 
                             B_chain  = B_chain[index_keep,], 
                             accept=accept/length(burnin:ttt), 
                             pscale=pscale, pcov = pcov, par_index=par_index)
            
            save(mcmc_out, file = paste0('Model_out/mcmc_out_',ind,'_', 'it', 
                                         ttt/chain_length_MASTER, '_samp', sampling_num, 
                                         '_', states_per_step, '_', steps_per_it,'.rda'))
            # Reset the chains
            chain = matrix( NA, chain_length_MASTER, length(par)) 
            B_chain = matrix( NA, chain_length_MASTER, length(y))
        }
    }
    # ---------------------------------------------------------------------------
    
    return(list( chain=chain, B_chain=B_chain, accept=accept/(steps-burnin), 
                 pscale=pscale, pcov = pcov, par_index=par_index))
}
