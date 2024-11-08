library(RcppArmadillo)
library(RcppDist)
library(Rcpp)

Rcpp::sourceCpp("likelihood_fnc.cpp")

dir = 'Model_out/'

trialNum = 2
sampNum = 2
itNum = 1
index_seeds = c(1:3)
long_chain = T

# -----------------------------------------------------------------------------
# Load the data ---------------------------------------------------------------
# -----------------------------------------------------------------------------
par_index = list()
par_index$vec_beta = 1:4
par_index$vec_alpha_tilde = 5:24
par_index$vec_sigma_upsilon = 25:424
par_index$vec_A = 425:444
par_index$vec_R = 445:460
par_index$vec_zeta = 461:484
par_index$vec_init = 485:488
par_index$omega_tilde = 489:576
par_index$vec_upsilon_omega = 577:664

sim_dat_num = 5

load('Data_sim/Dn_omega_names1.rda')
load('Data_sim/hr_map_names1.rda')
load(paste0('Data_sim/use_data1_', sim_dat_num, '.rda'))

Y = use_data[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(use_data[,'EID']))
eids = as.numeric(use_data[,'EID'])

x = use_data[,c('n_RBC_admin'), drop=F]
z = cbind(1, use_data[,c('RBC_ordered'), drop=F])

load(paste0('Data_sim/true_pars_', sim_dat_num, '.rda'))
true_par = true_pars

load(paste0('Data_sim/alpha_i_mat_', sim_dat_num, '.rda'))
load(paste0('Data_sim/omega_i_mat_', sim_dat_num, '.rda'))
load(paste0('Data_sim/Dn_omega_sim_', sim_dat_num, '.rda'))
load(paste0('Data_sim/bleed_indicator_sim_', sim_dat_num,'.rda'))

Dn_omega = Dn_omega_sim

b_chain = use_data[, "b_true"]

A = W = B_true = list()
for(i in EIDs){
    A[[i]] = alpha_i_mat[[which(EIDs == i)]]
    W[[i]] = omega_i_mat[[which(EIDs == i)]]
    B_true[[i]] = matrix(b_chain[eids == as.numeric(i)], ncol = 1)
}    

Dn_Xn_true = update_Dn_Xn_cpp( as.numeric(EIDs), B_true, Y, true_par, par_index, x, 5)
Dn_true = Dn_Xn_true[[1]]; names(Dn_true) = EIDs
Xn = Dn_Xn_true[[2]]

log_post_truth = log_post_cpp( as.numeric(EIDs), true_par, par_index, A, B_true,
                               Y, z, Dn_true, Xn, Dn_omega, W, 5)

# -----------------------------------------------------------------------------
# Load the MCMC chain and state sequence chains -------------------------------
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))
chain_list_b = vector(mode = "list", length = length(index_seeds))
ind = 0

for(seed in index_seeds){
    ind = ind + 1
    if(long_chain) {
        it_seq = 1:itNum
    } else {
        it_seq = itNum
    }
    
    for(it in it_seq) {
        
        file_name = paste0(dir,'mcmc_out_',toString(seed),'_', trialNum,'it', 
                           it, '_samp', sampNum, '_sim.rda') 
        
        load(file_name)
        print(paste0(ind, ": ", file_name))
        print("accept")
        print(mcmc_out_temp$accept)
        print("pscale")
        print(mcmc_out_temp$pscale)
        
        if(it == 1) {
            chain_list[[ind]] = mcmc_out_temp$chain
            chain_list_b[[ind]] = mcmc_out_temp$B_chain
        } else {
            chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out_temp$chain)
            chain_list_b[[ind]] = rbind(chain_list_b[[ind]], mcmc_out_temp$B_chain)
        }
        
        rm(mcmc_out_temp)
    }
}

EID_list = as.list(as.numeric(EIDs))
make_b = function(e, b_vec, eids) {
    return(matrix(b_vec[eids == e], ncol = 1))
}

keep_ind = seq(from = 5, to = nrow(chain_list[[1]]), by = 5)

log_post_chain = vector(mode = 'list', length = length(index_seeds))
for(i in 1:length(log_post_chain)) {
    log_post_chain[[i]] = rep(0, length(keep_ind))
}

for(i in 1:length(log_post_chain)) {
    for(j in 1:length(keep_ind)) {
        print(j)
        
        ind_j = keep_ind[j]
        
        par = chain_list[[i]][ind_j, ]
        b = chain_list_b[[i]][ind_j, ]
        
        B = lapply(EID_list, make_b, b_vec = b, eids = eids)
        
        Dn_Xn = update_Dn_Xn_cpp( as.numeric(EIDs), B, Y, par, par_index, x, 5)
        Dn = Dn_Xn[[1]]; names(Dn) = EIDs
        Xn = Dn_Xn[[2]]
        
        log_post_chain[[i]][j] = log_post_cpp( as.numeric(EIDs), par, par_index, A, B,
                                       Y, z, Dn, Xn, Dn_omega, W, 5)
    }
}

ylim_val = c(log_post_truth, do.call( c, log_post_chain))

plot( NULL, ylab=NA, main="log posterior", xlim=c(1,length(keep_ind)),
     ylim = c(min(ylim_val), max(ylim_val)),
     xlab = log_post_truth)

for(seed in 1:length(chain_list)) {
    lines( log_post_chain[[seed]], type='l', col=seed)
}
abline(h = log_post_truth)

