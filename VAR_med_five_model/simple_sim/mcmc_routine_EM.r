library(mvtnorm, quietly=T)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
library(expm, quietly = T)
sourceCpp("mcmc_fnc.cpp")

mcmc_routine_EM = function(par, par_index, B, y, ids, steps, burnin, ind){
    
    EIDs = unique(ids)
    interm_steps = 2000
    chain_interm = matrix( 0, interm_steps, length(par)) 
    chain_length_MASTER = 20
    steps = chain_length_MASTER * interm_steps
    
    # Number of cores over which to parallelize --------------------------------
    n_cores = 7
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
    chain = matrix( 0, chain_length_MASTER+1, length(par)) 
    B_chain = matrix( 0, chain_length_MASTER+1, length(y)) 
    accept = rep( 0, n_group)
    
    # Initialize B to the MLE based on initial values --------------------------
    B = viterbi_alg(as.numeric(EIDs), par, par_index, y, ids, n_cores)
    
    # Evaluate log-likelihood before MH step -----------------------------------
    log_target_prev = log_post_cpp(as.numeric(EIDs), par, par_index, B,
                                   y, ids, n_cores)
    if(!is.finite(log_target_prev)){
        print("Infinite log-posterior")
        print(paste0("value: ", log_target_prev))
        return(0)
    }
    
    # Metropolis-Hastings Algorithm --------------------------------------------
    B_chain[1, ] = do.call( 'c', B)
    chain[1, ] = par
    print("Current parameter values: "); print(par)
    
    chain_interm[1,] = par
    mcmc_start_t = Sys.time()
    for(ttt in 2:steps){
        
        if(ttt %% interm_steps == 0) {
            chain_ind = interm_steps
        } else {
            chain_ind = ttt - interm_steps * floor(ttt / interm_steps)
        }
        
        # Metropolis-Hastings updates ------------------------------------------
        for(j in 1:n_group) {
            
            ind_j = mpi[[j]]
            proposal = par
            
            proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                       sigma=pscale[[j]]*pcov[[j]])
            
            # Enforce mu1 > mu2
            if(sum(ind_j %in% par_index$mu) == length(ind_j)) {
                while(proposal[ind_j[1]] < proposal[ind_j[2]]) {
                    proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                               sigma=pscale[[j]]*pcov[[j]])
                }
            }
            
            # Evaluate proposed log-likelihood ---------------------------------
            log_target = log_post_cpp( as.numeric(EIDs), proposal, par_index, 
                                       B, y, ids, n_cores) 
            
            if(ttt < burnin){
                while(!is.finite(log_target)){
                    print('bad proposal')
                    proposal = par
                    proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                               sigma=pcov[[j]]*pscale[j])
                    
                    # Enforce mu1 > mu2
                    if(sum(ind_j %in% par_index$mu) == length(ind_j)) {
                        while(proposal[ind_j[1]] < proposal[ind_j[2]]) {
                            proposal[ind_j] = rmvnorm( n=1, mean=par[ind_j],
                                                       sigma=pscale[[j]]*pcov[[j]])
                        }
                    }
                    
                    log_target = log_post_cpp( as.numeric(EIDs), proposal, 
                                               par_index, B, y, ids, n_cores)
                }
            }
            
            if(!is.finite(log_target) | is.nan(log_target)) {
                # Ensuring that we do not have problems from C++
                print(paste0("bad proposal post burnin: ", log_target))
                log_target = -Inf
            }
            
            # Compute Metropolis ratio -----------------------------------------
            if( log_target - log_target_prev > log(runif(1,0,1)) ){
                log_target_prev = log_target
                par[ind_j] = proposal[ind_j]
                accept[j] = accept[j] +1
            }
            
            chain_interm[chain_ind,ind_j] = par[ind_j]
            
            # Proposal tuning scheme -------------------------------------------
            # During the burnin period, update the proposal covariance in each step
            # to capture the relationships within the parameters vectors for each
            # transition.  This helps with mixing.
            if(ttt < burnin){
                if(ttt == 100)  pscale[j] = 1
                
                if(100 <= ttt & ttt <= 2000){
                    temp_chain = chain_interm[1:ttt,ind_j]
                    pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
                    
                } else if(2000 < ttt){
                    temp_chain = chain_interm[(ttt-2000):ttt,ind_j]
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
            # ------------------------------------------------------------------
        }
        
        # Restart the acceptance ratio at burnin
        if(ttt == burnin) accept = rep( 0, n_group)
        # ----------------------------------------------------------------------
        
        if(ttt%%500==0)  print(paste0('--->',ttt, ' | ', log_target_prev))
        
        if(ttt > burnin & ttt%%interm_steps==0) {
            mcmc_end_t = Sys.time() - mcmc_start_t; print(mcmc_end_t)
            
            par_meds = apply(chain_interm, 2, median)
            par = par_meds
            B = mle_state_seq(as.numeric(EIDs), par, par_index, y, ids, n_cores)
            
            chain[ttt/interm_steps + 1, ] = par_meds
            B_chain[ttt/interm_steps + 1, ] = do.call( 'c', B)
            
            # Reset the chain
            chain_interm = matrix( NA, interm_steps, length(par)) 
            print("Current parameter values: "); print(par)
        }
    }
    # --------------------------------------------------------------------------
    
    return(list( chain=chain, B_chain=B_chain, accept=accept/(steps-burnin), 
                 pscale=pscale, pcov = pcov, par_index=par_index))
}
