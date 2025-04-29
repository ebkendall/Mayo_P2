args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])
simulation = as.logical(as.numeric(args[2]))

trialNum = 1
it_seq = 1:10
index_seeds = c(1:5)
states_per_step = 0
steps_per_it = 1

for(s in 1:length(index_seeds)) {
    seed_num = index_seeds[s]
    for(it in it_seq) {
        file_name = NULL
        if(simulation) {
            file_name = paste0('Model_out/mcmc_out_', trialNum,'_', seed_num,
                               'it', it, '_samp', sampling_num, '_', states_per_step,
                               '_', steps_per_it,'_sim.rda')
        } else {
            file_name = paste0('Model_out/mcmc_out_', trialNum,'_', seed_num,
                               'it', it, '_samp', sampling_num, '_', states_per_step,
                               '_', steps_per_it,'.rda')
        }
        
        if(file.exists(file_name)) {
            print(file_name)
            load(file_name)
            
            ind_keep = seq(2, nrow(mcmc_out_temp$chain), by = 2)
            
            mcmc_out_temp$chain    = mcmc_out_temp$chain[ind_keep, ]
            mcmc_out_temp$B_chain  = mcmc_out_temp$B_chain[ind_keep, ]
            mcmc_out_temp$hc_chain = mcmc_out_temp$hc_chain[ind_keep, ]
            mcmc_out_temp$hr_chain = mcmc_out_temp$hr_chain[ind_keep, ]
            mcmc_out_temp$bp_chain = mcmc_out_temp$bp_chain[ind_keep, ]
            mcmc_out_temp$la_chain = mcmc_out_temp$la_chain[ind_keep, ]
            
            for(a in 1:length(mcmc_out_temp$A_chain)) {
                mcmc_out_temp$A_chain[[a]] = mcmc_out_temp$A_chain[[a]][, ind_keep]
                if("W_chain" %in% names(mcmc_out_temp)) {
                    if(a == 1) print("omega RE exist")
                    mcmc_out_temp$W_chain[[a]] = mcmc_out_temp$W_chain[[a]][, ind_keep]
                }
            }
            
            save_file = NULL
            if(simulation) {
                save_file = paste0('Model_out/new/mcmc_out_', trialNum,'_', seed_num,
                                   'it', it, '_samp', sampling_num, '_', states_per_step,
                                   '_', steps_per_it,'_sim.rda')
            } else {
                save_file = paste0('Model_out/new/mcmc_out_', trialNum,'_', seed_num,
                                   'it', it, '_samp', sampling_num, '_', states_per_step,
                                   '_', steps_per_it,'.rda')
            }
            
            save(mcmc_out_temp, file = save_file)
            
            rm(mcmc_out_temp)
            
        } else {
            print(paste0("Missing! ", file_name))
        }
        
    }
}





# # Compile all results into a simple list ---------------------------------------
# 
# # Parameter estimates, pcov, pscale, state samples, imputation
# mcmc_results = vector(mode = 'list', length = length(index_seeds))
# 
# # Random effects
# a_w_chain_id = c(1, 86, 163, 237, 427)
# med_check_inds = c(2, 6, 9, 12, 29, 42, 46, 74)
# rand_eff_alpha = vector(mode = 'list', length = length(a_w_chain_id))
# # rand_eff_omega = vector(mode = 'list', length = length(a_w_chain_id))
# for(a in 1:length(a_w_chain_id)) {
#     rand_eff_alpha[[a]] = vector(mode = 'list', length = length(index_seeds))
#     # rand_eff_omega[[a]] = vector(mode = 'list', length = length(index_seeds))
# }
# 
# for(s in 1:length(index_seeds)) {
#     
#     seed_num = index_seeds[s]
#     
#     mcmc_results[[s]] = list()
#     mcmc_results[[s]]$pcov_pscale = list()
#     mcmc_results[[s]]$par = NULL
#     mcmc_results[[s]]$states = NULL
#     mcmc_results[[s]]$hemo = NULL
#     mcmc_results[[s]]$hr = NULL
#     mcmc_results[[s]]$map = NULL
#     mcmc_results[[s]]$lact = NULL
#     
#     pp_ind = 1
#     for(it in it_seq) {
#         
#         file_name = NULL
#         if(simulation) {
#             file_name = paste0('Model_out/mcmc_out_', trialNum,'_', seed_num,
#                                'it', it, '_samp', sampling_num, '_', states_per_step,
#                                '_', steps_per_it,'_sim.rda')
#         } else {
#             file_name = paste0('Model_out/mcmc_out_', trialNum,'_', seed_num,
#                                'it', it, '_samp', sampling_num, '_', states_per_step,
#                                '_', steps_per_it,'.rda')
#         }
#         
#         print(file_name)
#         load(file_name)
#         
#         if(pp_ind == 1) {
#             mcmc_results[[s]]$pcov_pscale[[1]] = mcmc_out_temp$pcov
#             mcmc_results[[s]]$pcov_pscale[[2]] = mcmc_out_temp$pscale
#         }
#         
#         ind_keep = seq(2, nrow(mcmc_out_temp$chain), by = 2)
#         
#         mcmc_results[[s]]$par    = rbind(mcmc_results[[s]]$par, mcmc_out_temp$chain[ind_keep, ])
#         mcmc_results[[s]]$states = rbind(mcmc_results[[s]]$states, mcmc_out_temp$B_chain[ind_keep, ])
#         mcmc_results[[s]]$hemo   = rbind(mcmc_results[[s]]$hemo, mcmc_out_temp$hc_chain[ind_keep, ])
#         mcmc_results[[s]]$hr     = rbind(mcmc_results[[s]]$hr, mcmc_out_temp$hr_chain[ind_keep, ])
#         mcmc_results[[s]]$map    = rbind(mcmc_results[[s]]$map, mcmc_out_temp$bp_chain[ind_keep, ])
#         mcmc_results[[s]]$lact   = rbind(mcmc_results[[s]]$lact, mcmc_out_temp$la_chain[ind_keep, ])
#         
#         for(a in 1:length(a_w_chain_id)) {
#             rand_eff_alpha[[a]][[s]] = cbind(rand_eff_alpha[[a]][[s]], mcmc_out_temp$A_chain[[a]][,ind_keep])
#             # rand_eff_omega[[a]][[s]] = cbind(rand_eff_omega[[a]][[s]], mcmc_out_temp$W_chain[[a]][,ind_keep])
#         }
#         
#         rm(mcmc_out_temp)
#         pp_ind = pp_ind + 1
#     }    
# }
# 
# save_file1 = NULL
# save_file2 = NULL
# save_file3 = NULL
# 
# if(simulation) {
#     save_file1 = paste0('Model_out/mcmc_results_', trialNum, '_samp', 
#                        sampling_num, '_it', max(it_seq), '_sim.rda')
#     save_file2 = paste0('Model_out/re_alpha_', trialNum, '_samp', 
#                         sampling_num, '_it', max(it_seq), '_sim.rda')
#     save_file3 = paste0('Model_out/re_omega_', trialNum, '_samp', 
#                         sampling_num, '_it', max(it_seq), '_sim.rda')
# } else {
#     save_file1 = paste0('Model_out/mcmc_results_', trialNum, '_samp', 
#                        sampling_num, '_it', max(it_seq), '.rda')
#     save_file2 = paste0('Model_out/re_alpha_', trialNum, '_samp', 
#                         sampling_num, '_it', max(it_seq), '.rda')
#     save_file3 = paste0('Model_out/re_omega_', trialNum, '_samp', 
#                         sampling_num, '_it', max(it_seq), '.rda')
# }
# 
# save(mcmc_results, file = save_file1)
# save(rand_eff_alpha, file = save_file2)
# # save(rand_eff_omega, file = save_file3)



