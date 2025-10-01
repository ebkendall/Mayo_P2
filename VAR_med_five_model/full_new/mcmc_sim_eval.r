index_seeds = c(1:25)

trialNum = 1
S = 5

it_num = 2
start_ind = 1

# Mode of the state sequences -------------------------------------------------
Mode <- function(x) {
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
}

# Load the model output -------------------------------------------------------
state_results = list()
seed_list = NULL

ind = 0
for(seed in index_seeds){
    
    it_seq = 1:(it_num - start_ind + 1)
    
    B_chain = NULL
    
    check_name = paste0('Model_out/mcmc_out_',trialNum, '_', seed, 'it', it_num,'_sim.rda')    
    
    if(file.exists(check_name)) {
        
        for(it in it_seq) {
            
            load(paste0('Model_out/mcmc_out_',trialNum, '_', seed, 
                        'it', it + (start_ind - 1),'_sim.rda'))
            
            print(paste0(seed, ": ", it))
            par_index = mcmc_out$par_index
            
            if(it == 1) {
                ind = ind + 1
                B_chain = mcmc_out$B_chain[250:500, ]
            } else {
                B_chain = rbind(B_chain, mcmc_out$B_chain)
            }
            
            rm(mcmc_out)
        }
        
        state_results[[ind]] = list()
        state_results[[ind]][[1]] = matrix(nrow = S+1, ncol = ncol(B_chain))
        for(jj in 1:S) {
            state_results[[ind]][[1]][jj, ] = apply(B_chain, 2, function(x,jj){sum(x == jj)}, jj)
        }
        state_results[[ind]][[1]][S+1, ] = apply(B_chain, 2, Mode)
        
        seed_list = c(seed_list, seed)
        
        state_2_cumulative = colMeans(B_chain == 2)
        
        state_results[[ind]][[2]] = state_2_cumulative
        
    } else {
        print(paste0("Missing! ", check_name))
    }
    
    rm(B_chain)
}

# save(state_results, file = paste0("Model_out/state_results_sim_", it_num, ".rda"))
 
# Computing ROC curve information ---------------------------------------------
mode_correct = rep(NA, length(state_results))
c = seq(0, 1, by = 0.001)

chosen_c = chosen_c_alt = NULL
chosen_c_ppv = chosen_c_ppv_alt = NULL
computed_AUC = computed_AUC_alt = rep(NA, length(seed_list))

pdf(paste0('Plots/ROC_og_it', it_num, '.pdf'))
par(mfrow=c(3, 2))
for(ii in 1:length(seed_list)) {
    
    seed = seed_list[ii]
    
    print(seed)

    load(paste0('Data/sim_data_', seed, '.rda'))
    EIDs = unique(data_format[,"EID"])
    true_state = data_format[,'b_true']

    # Quick summary
    mode_correct[ii] = mean(state_results[[ii]][[1]][S+1, ] == true_state)

    # State 2 focus (OG method) --------------------
    state_2_binary = as.numeric(true_state == 2)
    state_1345_binary = as.numeric(true_state != 2)
    
    cumulative_post_prob = state_results[[ii]][[2]]
    
    success_mat = matrix(nrow=length(c), ncol = length(state_2_binary))
    for(i in 1:length(c)) {
        success_mat[i,] = as.numeric(cumulative_post_prob >= c[i])
    }
    
    true_positives = which(state_2_binary == 1)
    true_negatives = which(state_1345_binary == 1)
    
    p = length(true_positives)  # Number of true "positives"
    n = length(true_negatives)  # Number of true "negatives"
    
    c_results = data.frame("c" = c, "true_pos"  = rep(NA, length(c)),
                           "false_pos" = rep(NA, length(c)),
                           "true_neg"  = rep(NA, length(c)),
                           "false_neg" = rep(NA, length(c)),
                           "poss_pred" = rep(NA, length(c)))
    
    # State 2 focus (Alt method) --------------------
    state_2_binary_alt = rep(0, length(true_state))
    for(i in 1:length(EIDs)) {
        eid = EIDs[i]
        eid_ind = which(data_format[,"EID"] == eid)
        sub_df_true_state = data_format[eid_ind, "b_true"]
        
        s2_indices = which(sub_df_true_state == 2)
        if(length(s2_indices) > 0) {
            s2_indices_new = sort(unique(c(s2_indices - 1, s2_indices, s2_indices+1)))
            s2_indices_new = s2_indices_new[s2_indices_new > 0 & s2_indices_new <= length(eid_ind)]
            
            s2_indicator = rep(0, length(eid_ind))
            s2_indicator[s2_indices_new] = 1
            
            state_2_binary_alt[eid_ind] = s2_indicator    
        }
    }
    
    true_positives_alt = which(state_2_binary_alt == 1)
    true_negatives_alt = which(state_2_binary_alt == 0)
    
    p_alt = length(true_positives_alt)  # Number of true "positives"
    n_alt = length(true_negatives_alt)  # Number of true "negatives"
    
    c_results_alt = data.frame("c" = c, "true_pos"  = rep(NA, length(c)),
                               "false_pos" = rep(NA, length(c)),
                               "true_neg"  = rep(NA, length(c)),
                               "false_neg" = rep(NA, length(c)),
                               "poss_pred" = rep(NA, length(c)))
    
    for(i in 1:nrow(success_mat)) {
        
        test_positives = which(success_mat[i,] == 1)
        test_negatives = which(success_mat[i,] == 0)
        
        # original AUC approach -------------------------
        tpr = sum(test_positives %in% true_positives) / p
        tnr = sum(test_negatives %in% true_negatives) / n
        
        fnr = 1 - tpr
        fpr = 1 - tnr
        
        tp = sum(test_positives %in% true_positives)
        fp = sum(test_positives %in% true_negatives)
        
        ppv = tp / (tp + fp)
        
        c_results[i,] = c(c[i], tpr, fpr, tnr, fnr, ppv)
        
        # alternate AUC approach -------------------------
        tpr_alt = sum(test_positives %in% true_positives_alt) / p_alt
        tnr_alt = sum(test_negatives %in% true_negatives_alt) / n_alt
        
        fnr_alt = 1 - tpr_alt
        fpr_alt = 1 - tnr_alt
        
        tp_alt = sum(test_positives %in% true_positives_alt)
        fp_alt = sum(test_positives %in% true_negatives_alt)
        
        ppv_alt = tp_alt / (tp_alt + fp_alt)
        
        c_results_alt[i,] = c(c[i], tpr_alt, fpr_alt, tnr_alt, fnr_alt, ppv_alt)
        
    }
    
    # original AUC approach -------------------------
    temp = c_results[order(c_results$false_pos), ]
    c_results = temp
    
    min_diff = c_results$false_pos - c_results$true_pos
    poss_c = c_results$c[which(min_diff == min(min_diff))]
    
    max_ppv = c_results$poss_pred
    poss_c_ppv = c_results$c[which(max_ppv == max(max_ppv))]
    
    chosen_c = c(chosen_c, poss_c)
    chosen_c_ppv = c(chosen_c_ppv, poss_c_ppv)
    
    # AUC (estimated using Trapezoid rule)
    area = 0
    for(k in 2:nrow(c_results)) {
        area = area + 0.5 * (c_results$true_pos[k] + c_results$true_pos[k-1]) *
            (c_results$false_pos[k] - c_results$false_pos[k-1])
    }
    
    computed_AUC[ii] = area
    
    grid = seq(0,1,length = nrow(c_results))
    s1 <- smooth.spline(grid, c_results$false_pos, all.knots = T, penalty = 0)
    s2 <- smooth.spline(grid, c_results$true_pos, all.knots = T, penalty = 0)
    
    x2 <- predict(s1, grid)[[2]]
    y2 <- predict(s2, grid)[[2]]
    
    plot(c_results$false_pos, c_results$true_pos, xlim = c(0,1), ylim = c(0,1),
         xlab = "FPR", ylab = "TPR",
         main = paste0("OG: Seed = ", seed, ", AUC = ", round(area, digits = 4)))
    points(c_results$false_pos[c_results$c %in% poss_c],
           c_results$true_pos[c_results$c %in% poss_c],
           col = 'red', cex = 2)
    points(c_results$false_pos[c_results$c %in% poss_c_ppv],
           c_results$true_pos[c_results$c %in% poss_c_ppv],
           col = 'purple', cex = 2)
    abline(a = 0, b = 1, col = 'red', lty = 2)
    lines(x2, y2, col= 'green')
    
    # alternate AUC approach -------------------------
    temp = c_results_alt[order(c_results_alt$false_pos), ]
    c_results_alt = temp
    
    min_diff = c_results_alt$false_pos - c_results_alt$true_pos
    poss_c = c_results_alt$c[which(min_diff == min(min_diff))]
    
    max_ppv = c_results_alt$poss_pred
    poss_c_ppv = c_results_alt$c[which(max_ppv == max(max_ppv))]
    
    chosen_c_alt = c(chosen_c_alt, poss_c)
    chosen_c_ppv_alt = c(chosen_c_ppv_alt, poss_c_ppv)
    
    # AUC (estimated using Trapezoid rule)
    area = 0
    for(k in 2:nrow(c_results_alt)) {
        area = area + 0.5 * (c_results_alt$true_pos[k] + c_results_alt$true_pos[k-1]) *
            (c_results_alt$false_pos[k] - c_results_alt$false_pos[k-1])
    }
    
    computed_AUC_alt[ii] = area
    
    grid = seq(0,1,length = nrow(c_results_alt))
    s1 <- smooth.spline(grid, c_results_alt$false_pos, all.knots = T, penalty = 0)
    s2 <- smooth.spline(grid, c_results_alt$true_pos, all.knots = T, penalty = 0)
    
    x2 <- predict(s1, grid)[[2]]
    y2 <- predict(s2, grid)[[2]]
    
    plot(c_results_alt$false_pos, c_results_alt$true_pos, xlim = c(0,1), ylim = c(0,1),
         xlab = "FPR", ylab = "TPR",
         main = paste0("ALT: Seed = ", seed, ", AUC = ", round(area, digits = 4)))
    points(c_results_alt$false_pos[c_results_alt$c %in% poss_c],
           c_results_alt$true_pos[c_results_alt$c %in% poss_c],
           col = 'red', cex = 2)
    points(c_results_alt$false_pos[c_results_alt$c %in% poss_c_ppv],
           c_results_alt$true_pos[c_results_alt$c %in% poss_c_ppv],
           col = 'purple', cex = 2)
    abline(a = 0, b = 1, col = 'red', lty = 2)
    lines(x2, y2, col= 'green')
    
}
dev.off()

# save(chosen_c, file = paste0("Model_out/chosen_c_sim_", it_num, ".rda"))
# save(chosen_c_ppv, file = paste0("Model_out/chosen_c_ppv_sim_", it_num, ".rda"))
# save(computed_AUC, file = paste0("Model_out/AUC_sim_", it_num, ".rda"))

print("Summary of state identification based on posterior mode: ")
print(summary(mode_correct))

print("ORIGINAL METHOD")

print("Choice for c (FPR - TPR)")
print(summary(chosen_c))
indiv_select_c = median(chosen_c)

print("Choice for c (PPV)")
print(summary(chosen_c_ppv))

print("Average AUC")
print(summary(computed_AUC))

sens_and_spec_S2 = matrix(nrow = length(seed_list), ncol = 2)
colnames(sens_and_spec_S2) = c("sens", "spec")
for(ii in 1:length(seed_list)) {
    
    seed = seed_list[ii]
    
    load(paste0('Data/sim_data_', seed, '.rda'))
    true_state = data_format[,'b_true']

    post_prob_S2 = state_results[[ii]][[1]][2, ] / colSums(state_results[[ii]][[1]][1:S,])
    S2_identification = as.numeric(post_prob_S2 >= indiv_select_c)

    # Sensitivity of state 2: Pr(predict S2 | true S2)
    true_S2_ind = S2_identification[true_state == 2]
    sens_S2 = mean(true_S2_ind == 1)

    # Specificity of state 2: Pr(predict not S2 | true not S2)
    true_not_S2_ind = S2_identification[true_state != 2]
    spec_S2 =  mean(true_not_S2_ind == 0)

    sens_and_spec_S2[ii, ] = c(sens_S2, spec_S2)
}

print(paste0("Sensitivity of S2 using c = ", indiv_select_c))
print(summary(sens_and_spec_S2[,"sens"]))

print(paste0("Specificity of S2 using c = ", indiv_select_c))
print(summary(sens_and_spec_S2[,"spec"]))

print("ALTERNATE METHOD")

print("Choice for c (FPR - TPR)")
print(summary(chosen_c_alt))
indiv_select_c_alt = median(chosen_c_alt)

print("Choice for c (PPV)")
print(summary(chosen_c_ppv_alt))

print("Average AUC")
print(summary(computed_AUC_alt))

sens_and_spec_S2 = matrix(nrow = length(seed_list), ncol = 2)
colnames(sens_and_spec_S2) = c("sens", "spec")
for(ii in 1:length(seed_list)) {
    
    seed = seed_list[ii]
    
    load(paste0('Data/sim_data_', seed, '.rda'))
    true_state = data_format[,'b_true']
    
    post_prob_S2 = state_results[[ii]][[1]][2, ] / colSums(state_results[[ii]][[1]][1:S,])
    S2_identification = as.numeric(post_prob_S2 >= indiv_select_c_alt)
    
    # Sensitivity of state 2: Pr(predict S2 | true S2)
    true_S2_ind = S2_identification[true_state == 2]
    sens_S2 = mean(true_S2_ind == 1)
    
    # Specificity of state 2: Pr(predict not S2 | true not S2)
    true_not_S2_ind = S2_identification[true_state != 2]
    spec_S2 =  mean(true_not_S2_ind == 0)
    
    sens_and_spec_S2[ii, ] = c(sens_S2, spec_S2)
}

print(paste0("Sensitivity of S2 using c = ", indiv_select_c_alt))
print(summary(sens_and_spec_S2[,"sens"]))

print(paste0("Specificity of S2 using c = ", indiv_select_c_alt))
print(summary(sens_and_spec_S2[,"spec"]))


