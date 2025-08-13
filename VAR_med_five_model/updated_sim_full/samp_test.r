library(mvtnorm, quietly=T)
library(LaplacesDemon, quietly=T)
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
sourceCpp("samp_test.cpp")

# # Needed for OpenMP C++ parallel
# Sys.setenv("PKG_CXXFLAGS" = "-fopenmp")
# Sys.setenv("PKG_LIBS" = "-fopenmp")

set.seed(2025)
load('Data/sim_data_1.rda')
load('Data/alpha_i_mat_1.rda')
load('Data/omega_i_mat_1.rda')
load('Data/bleed_indicator_sim_1.rda')
load('Data/Dn_omega_sim.rda')
EIDs_all = unique(data_format[,'EID'])

# Take 10 simulated subjects ---------------------------------------------------
EIDs = sort(sample(EIDs_all, size = 10, replace = F))

index_sub = (data_format[,"EID"] %in% EIDs)

Y = data_format[index_sub, c('EID', 'hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
x = data_format[index_sub, c('n_RBC_admin'), drop=F]
z = cbind(1, data_format[index_sub, c('RBC_ordered'), drop=F])

bleed_indicator = bleed_indicator[index_sub]
Dn_omega_sim = Dn_omega_sim[which(EIDs_all %in% EIDs)]
Dn_omega = Dn_omega_sim; rm(Dn_omega_sim)
alpha_i_mat = alpha_i_mat[which(EIDs_all %in% EIDs)]
omega_i_mat = omega_i_mat[which(EIDs_all %in% EIDs)]

true_b_chain = data_format[index_sub, "b_true"]
b_chain = rep(1, sum(index_sub))

A = list()
W = list()
B = list()

for(ii in 1:length(EIDs)){
    i = EIDs[ii]
    
    A[[ii]] = alpha_i_mat[[ii]][-c(1,6,11,16),,drop=F] # Remove baseline
    W[[ii]] = omega_i_mat[[ii]]
    B[[ii]] = matrix(b_chain[Y[,"EID"] == i], ncol = 1)
}

Xn = initialize_Xn(EIDs, Y, x)
Dn_alpha = initialize_Dn(EIDs, B)

par_index = list()
par_index$beta = 1:4
par_index$alpha_tilde = 5:20
par_index$upsilon = 21:276
par_index$A = 277:280
par_index$R = 281:296
par_index$zeta = 297:320
par_index$init = 321:324
par_index$omega_tilde = 325:408
par_index$eta_omega = 409:492
par_index$G = 493:508

par = rep(0, max(do.call('c', par_index)))

par[par_index$beta] = c(0.25, -2, 2, -0.25) 
par[par_index$alpha_tilde] = c( -5,   5, -2,  2,
                                10, -10,  2, -2,
                               -10,  10,  2, -2,
                                 5,  -5, -2,  2)
par[par_index$upsilon] = c(diag(c(4, 4, 1, 1, 
                                  4, 4, 1, 1, 
                                  4, 4, 1, 1, 
                                  4, 4, 1, 1)))
par[par_index$A] = rep(0, 4)
par[par_index$R] = c(diag(c(9, 9, 9, 9)))
#    transitions:          1->2,         1->4,         2->3,         2->4, 
#                          3->1,         3->2,         3->4,         4->2, 
#                          4->5,         5->1,         5->2,         5->4
par[par_index$zeta] = c(-3.7405, 2.5, -4.2152,   1, -2.6473,-0.5, -2.1475, -0.2, 
                        -3.4459,  -1, -2.9404,   1, -3.2151,   1, -3.1778,  1.5, 
                        -2.0523,   0, -3.4459,-0.2, -3.2404, 2.5, -3.2151,    1)

init_dist = c(160,  82,  38, 104, 116)
init_dist = init_dist / sum(init_dist)
par[par_index$init[1]] = log(init_dist[2] / (1 - sum(init_dist[2:5])))
par[par_index$init[2]] = log(init_dist[3] / (1 - sum(init_dist[2:5])))
par[par_index$init[3]] = log(init_dist[4] / (1 - sum(init_dist[2:5])))
par[par_index$init[4]] = log(init_dist[5] / (1 - sum(init_dist[2:5])))

par[par_index$omega_tilde]= 2 * c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                                  -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                                  -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                                  -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
par[par_index$eta_omega] = rep(1, length(par_index$eta_omega))
par[par_index$G] = c(diag(c(9, 9, 9, 9)))


for(samp_ind in 1:4) {
    
    if(samp_ind == 1) {
        sps = 2
    } else if(samp_ind == 2) {
        sps = 2
    } else if(samp_ind == 3) {
        sps = 2
    } else {
        sps = 2       
    }
    
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
    
    initialize_cpp(adj_mat, adj_mat_sub, sps)
    
    B_chain = matrix(NA, 1000, nrow(Y)) 
    B_chain[1, ] = do.call( 'c', B)
    
    for(ttt in 2:steps) {
        
        # Random sample update
        if(sampling_num == 1) {
            B_Dn = mh_up(EIDs, par, par_index, A, B, Y, z, Xn, Dn_omega, W, 
                         bleed_indicator, n_cores, states_per_step)
            B = B_Dn[[1]]
            Dn = B_Dn[[2]]
        }
        
        # Almost-Gibbs update
        else if(sampling_num == 2) {
            B_Dn = almost_gibbs_up(EIDs, par, par_index, A, B, Y, z, Xn, 
                                   Dn_omega, W, bleed_indicator, n_cores, 
                                   states_per_step)
            B = B_Dn[[1]]
            Dn = B_Dn[[2]]
        } 
        
        # Gibbs update
        else if(sampling_num == 3) {
            B_Dn = gibbs_up(EIDs, par, par_index, A, B, Y, z, Xn, Dn_omega,
                            W, bleed_indicator, n_cores, states_per_step)
            B = B_Dn[[1]]
            Dn = B_Dn[[2]]
        }
        
        # Full seq MH update
        else if(sampling_num == 4) {
            B_Dn = mh_up_all(EIDs, par, par_index, A, B, Y, z, Xn, Dn_omega,
                             W, bleed_indicator, n_cores)
            B = B_Dn[[1]]
            Dn = B_Dn[[2]]
        }
        
        # Almost-Gibbs efficient
        else if(sampling_num == 5) {
            
            if(states_per_step == 0) {
                sps = sample(x = 20:50, size = 1, replace = T)
            } else {
                sps = states_per_step
            }
            
            B_Dn = almost_gibbs_fast_b(EIDs, par, par_index, A, B, Y, z, Xn, 
                                       Dn_omega, W, bleed_indicator,n_cores,
                                       sps)
            B = B_Dn[[1]] 
            Dn = B_Dn[[2]]
        }
    }
}


