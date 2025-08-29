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

Mode <- function(x) {
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
}

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

for(ii in 1:length(EIDs)){
    A[[ii]] = alpha_i_mat[[ii]][-c(1,6,11,16),,drop=F] # Remove baseline
    W[[ii]] = omega_i_mat[[ii]]
}

Xn = initialize_Xn(EIDs, Y, x)

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

# Run time experiment ----------------------------------------------------------
compute_num_combos = cbind(rep(1:3, each = 5), rep(1:5, 3))
colnames(compute_num_combos) = c("samp_ind", "s")
compute_num_combos = rbind(compute_num_combos, c(4, 1))

args = commandArgs(TRUE)
compute_num = as.numeric(args[1])

samp_ind = compute_num_combos[compute_num, "samp_ind"]
s = compute_num_combos[compute_num, "s"]

compute_times = list()
steps = 100
n_cores = 1

# for(samp_ind in 1:4) {
    
    if(samp_ind < 4) {
        sps = c(2,4,6,8,10)
    } else {
        sps = 2
    }
    
    # compute_times[[samp_ind]] = list()
    
    # for(s in 1:length(sps)) {
        
        initialize_cpp(adj_mat, adj_mat_sub, sps[s])
        get_dimension()
        
        B = list()
        for(ii in 1:length(EIDs)){
            i = EIDs[ii]
            B[[ii]] = matrix(b_chain[Y[,"EID"] == i], ncol = 1)
        }
        
        Dn_alpha = initialize_Dn(EIDs, B)
        
        B_chain = matrix(NA, steps, nrow(Y)) 
        B_chain[1, ] = do.call( 'c', B)
        
        interm_compute_time = matrix(nrow = steps, ncol = 2)
        colnames(interm_compute_time) = c('time', 'accuracy')
        interm_compute_time[1,] = c(0, mean(B_chain[1,] == true_b_chain))
        
        # compute_times[[samp_ind]][[s]] = matrix(nrow = steps, ncol = 2)
        # colnames(compute_times[[samp_ind]][[s]]) = c('time', 'accuracy')
        # compute_times[[samp_ind]][[s]][1,] = c(0, mean(B_chain[1,] == true_b_chain))
        
        total_time_start = Sys.time()
        for(ttt in 2:steps) {
            
            ttt_start_t = Sys.time()
            
            # Coin Flip
            if(samp_ind == 1) {
                B_Dn = state_coin_flip(EIDs, par, par_index, B, A, W, Y, z, 
                                       Dn_omega, Xn, bleed_indicator, 
                                       sps[s], n_cores)
                B = B_Dn[[1]]
                Dn_alpha = B_Dn[[2]]
            }
            
            # Almost-Gibbs update
            else if(samp_ind == 2) {
                B_Dn = state_almost_gibbs(EIDs, par, par_index, B, A, W, Y, z, 
                                          Dn_omega, Xn, bleed_indicator, 
                                          sps[s], n_cores)
                B = B_Dn[[1]]
                Dn_alpha = B_Dn[[2]]
            } 
            
            # Gibbs update
            else if(samp_ind == 3) {
                B_Dn = state_gibbs(EIDs, par, par_index, B, A, W, Y, z, 
                                   Dn_omega, Xn, bleed_indicator, 
                                   sps[s], n_cores)
                B = B_Dn[[1]]
                Dn_alpha = B_Dn[[2]]
            }
            
            # Our Sampler
            else {
                state_per_step = sample(x = 2:50, size = 1, replace = T)
                B_Dn = state_sampler(EIDs, par, par_index, B, A, W, Y, z, 
                                     Dn_omega, Xn, bleed_indicator, 
                                     state_per_step, n_cores)
                B = B_Dn[[1]]
                Dn_alpha = B_Dn[[2]]
            }
            
            B_chain[ttt, ] = do.call( 'c', B)
            mode_chain = apply(B_chain[1:ttt, ], 2, Mode)
            
            ttt_end_t = Sys.time()
            
            ttt_elapsed = as.numeric(difftime(ttt_end_t, ttt_start_t, units = "secs"))
            ttt_accuracy = mean(mode_chain == true_b_chain)
            
            interm_compute_time[ttt,] = c(ttt_elapsed, ttt_accuracy)
            # compute_times[[samp_ind]][[s]][ttt,] = c(ttt_elapsed, mean(mode_chain == true_b_chain))
            
            cat('--> samp ind ', samp_ind, ', sps ', sps[s], ', ', ttt, '\n')
            cat("Elapsed time:", ttt_elapsed, "seconds\n")
            cat("Accuracy:", ttt_accuracy, "\n")
        }
        total_time_end = Sys.time()
        total_time_elapsed = as.numeric(difftime(total_time_end, total_time_start, units = "secs"))
        cat("Total elapsed time:", total_time_elapsed, "seconds\n")
    # }
# }
        
        
save(interm_compute_time, file = paste0('Model_out/int_comp_', samp_ind, '_', s, '.rda'))
# save(compute_times, file = 'Model_out/compute_times.rda')

compute_times = list()
for(i in 1:4) {
    compute_times[[i]] = list()
    for(s in 1:5) {
        file_name = paste0('Model_out/int_comp_', i, '_', s, '.rda')
        if(file.exists(file_name)) {
            load(file_name)
            compute_times[[i]][[s]] = interm_compute_time
        }
    }
}

compute_times_10 = c(1832.113, 2948.188, 3024.797, 2986.866, 3023.945, 3057.364,
                     2996.636, 2984.072, 2989.379, 3066.625, 3009.258, 3007.343, 
                     3044.969, 3107.056)
compute_times_14 = c(1724.982, 2850.964, 2891.257, 2853.593, 2845.207, 2803.919,
                     2788.594, 2795.35,  2791.955, 2798.976, 2790.874, 2787.029,
                     2858.222, 2926.53,  2755.723)
compute_times_15 = c(20454.7)

pdf("Plots/compute_times_update.pdf")
par(mfrow = c(3, 2))
plot_names = c("Coin Flip", "Almost-Gibbs", "Gibbs", "Our Sampler")
for(i in 1:4) {
    for(j in 1:length(compute_times[[i]])) {
        plot(compute_times[[i]][[j]][,2], main = paste0(plot_names[i], ", p = ", sps[j]),
             ylim = c(0,1), ylab = "Percent Correct",
             xlab = paste0("Median Time = ", round(median(compute_times[[i]][[j]][,1]), 4),
                           ", Accuracy = ", round(compute_times[[i]][[j]][100,2], 4)))
    }
}
dev.off()
