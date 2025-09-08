library(foreach, quietly=T)
library(doParallel, quietly=T)

index_seeds = c(1:20)

trialNum = 1
S = 5

it_num = 2
start_ind = 1

window_length = 0:11

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
    
    B_chain   = NULL
    
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
        
        # Testing different cumulative measures of posterior probabilities -----
        load(paste0('Data/sim_data_', seed, '.rda'))
        EIDs = unique(data_format[,"EID"])
        
        cores=detectCores()
        cl <- makeCluster(cores[1]-1) 
        registerDoParallel(cl)
        
        start_t = Sys.time()
        state_2_cumulative = foreach(www=1:length(window_length), .inorder = T) %dopar% {
            win_length = window_length[www]
            
            cumulative_post_prob = rep(NA, nrow(data_format))
            cumulative_post_prob_2 = rep(NA, nrow(data_format))
            ttt = 1
            for(i in EIDs) {
                print(which(EIDs == i))
                indices_i = which(data_format[,'EID']==i)
                for(w in 1:length(indices_i)) {
                    start_index = indices_i[1]
                    end_index = indices_i[w]
                    if(w - win_length > 0) start_index = indices_i[w - win_length]
                    
                    y_or_n_2 = apply(B_chain[, start_index:end_index, drop=F],
                                     1, function(x) (2 %in% x))
                    prob_2 = mean(y_or_n_2)
                    
                    cumulative_post_prob[ttt] = prob_2
                    
                    alt_y_or_n_2 = mean(c(B_chain[, start_index:end_index]) == 2)
                    cumulative_post_prob_2[ttt] = alt_y_or_n_2
                    
                    ttt = ttt + 1
                }
            }
            
            both_cumulative = rbind(cumulative_post_prob, cumulative_post_prob_2)
            
            return(both_cumulative)
        }
        end_t = Sys.time(); 
        ttt_elapsed = as.numeric(difftime(end_t, start_t, units = "secs"))
        print(ttt_elapsed)
        
        stopCluster(cl)
        
        state_results[[ind]][[2]] = state_2_cumulative[[1]][1,]
        state_results[[ind]][[3]] = state_2_cumulative[[1]][2,]
        
        for(k in 2:length(window_length)) {
            state_results[[ind]][[2]] = rbind(state_results[[ind]][[2]], state_2_cumulative[[k]][1,]) 
            state_results[[ind]][[3]] = rbind(state_results[[ind]][[3]], state_2_cumulative[[k]][2,])
        }
        
        rm(data_format)
        
    } else {
        print(paste0("Missing! ", check_name))
    }
    
    rm(B_chain)
}

save(state_results, file = paste0("Model_out/state_results_sim_", it_num, ".rda"))
 
# Computing ROC curve information ---------------------------------------------
mode_correct = rep(NA, length(state_results))
c = seq(0, 1, by = 0.001)
chosen_c = list()
seed_list = index_seeds
computed_AUC = matrix(nrow = length(seed_list), ncol = length(window_length))

pdf(paste0('Plots/ROC_it', it_num, '.pdf'))
par(mfrow=c(3, 2))
for(seed in seed_list) {
    print(seed)

    load(paste0('Data/sim_data_', seed, '.rda'))
    EIDs = unique(data_format[,"EID"])
    true_state = data_format[,'b_true']

    # Quick summary
    mode_correct[seed] = mean(state_results[[seed]][[1]][S+1, ] == true_state)

    # State 2 focus
    state_2_binary = as.numeric(true_state == 2)
    state_1345_binary = as.numeric(true_state != 2)

    chosen_c[[seed]] = list()

    for(www in 1:nrow(state_results[[seed]][[2]])) {

        print(www)

        cumulative_post_prob = state_results[[seed]][[2]][www, ]

        success_mat = matrix(nrow=length(c), ncol = length(state_2_binary))
        for(i in 1:length(c)) {
            success_mat[i,] = as.numeric(cumulative_post_prob >= c[i])
        }

        true_positives = which(state_2_binary == 1)
        true_negatives = which(state_1345_binary == 1)

        c_results = data.frame("c" = c, "true_pos"  = rep(NA, length(c)),
                               "false_pos" = rep(NA, length(c)),
                               "true_neg"  = rep(NA, length(c)),
                               "false_neg" = rep(NA, length(c)))

        # Calculating the sensitivity and specificity information
        p = length(true_positives)  # Number of true "positives"
        n = length(true_negatives)  # Number of true "negatives"

        for(i in 1:nrow(success_mat)) {
            test_positives = which(success_mat[i,] == 1)
            test_negatives = which(success_mat[i,] == 0)

            tpr = sum(test_positives %in% true_positives) / p
            tnr = sum(test_negatives %in% true_negatives) / n

            fnr = 1 - tpr
            fpr = 1 - tnr

            c_results[i,] = c(c[i], tpr, fpr, tnr, fnr)
        }

        temp = c_results[order(c_results$false_pos), ]
        c_results = temp

        min_diff = c_results$false_pos - c_results$true_pos
        poss_c = c_results$c[which(min_diff == min(min_diff))]

        chosen_c[[seed]][[www]] = poss_c

        # slope_mat = cbind(diff(c_results$true_pos), diff(c_results$false_pos))
        # slope_mat = cbind(c[1:nrow(slope_mat)], slope_mat)
        # colnames(slope_mat) = c("c", "delta_y", "delta_x")
        #
        # slope_mat = slope_mat[slope_mat[,"delta_x"] != 0, ]
        # slopes = slope_mat[,"delta_y"] / slope_mat[,"delta_x"]
        # slope_mat = cbind(slope_mat, slopes)
        #
        # approx_1_slope = slope_mat[order(abs(slope_mat[,"slopes"] - 1))[1], "slopes"]
        # c_slope_1 = slope_mat[slope_mat[,"slopes"] == approx_1_slope, "c"]
        # chosen_c = rbind(chosen_c, c_results[c_results$c %in% c_slope_1, ])

        # AUC (estimated using Trapezoid rule)
        area = 0
        for(k in 2:nrow(c_results)) {
            area = area + 0.5 * (c_results$true_pos[k] + c_results$true_pos[k-1]) *
                (c_results$false_pos[k] - c_results$false_pos[k-1])
        }

        computed_AUC[seed, www] = area

        grid = seq(0,1,length = nrow(c_results))
        s1 <- smooth.spline(grid, c_results$false_pos, all.knots = T, penalty = 0)
        s2 <- smooth.spline(grid, c_results$true_pos, all.knots = T, penalty = 0)

        x2 <- predict(s1, grid)[[2]]
        y2 <- predict(s2, grid)[[2]]

        plot(c_results$false_pos, c_results$true_pos, xlim = c(0,1), ylim = c(0,1),
             xlab = "FPR", ylab = "TPR",
             main = paste0("Seed = ", seed, ", AUC = ", round(area, digits = 4), ", w = ", window_length[www]))
        points(c_results$false_pos[c_results$c %in% chosen_c[[seed]][[www]]],
               c_results$true_pos[c_results$c %in% chosen_c[[seed]][[www]]],
               col = 'red', cex = 2)
        abline(a = 0, b = 1, col = 'red', lty = 2)
        lines(x2, y2, col= 'green')
    }
}
dev.off()

save(chosen_c, file = paste0("Model_out/chosen_c_sim_", it_num, ".rda"))
save(computed_AUC, file = paste0("Model_out/AUC_sim_", it_num, ".rda"))

print("Summary of state identification based on posterior mode: ")
print(summary(mode_correct))

print("Summary of the the choice for c")
all_c = NULL
for(i in 1:length(chosen_c)) {
    all_c = c(all_c, chosen_c[[i]][[1]])
}
print(summary(all_c))

print("Average AUC")
all_AUC = NULL
for(i in 1:length(chosen_c)) {
    all_AUC = c(all_AUC, computed_AUC[i,1])
}
print(summary(all_AUC))

indiv_select_c = mean(all_c)

sens_and_spec_S2 = matrix(nrow = length(seed_list), ncol = 2)
colnames(sens_and_spec_S2) = c("sens", "spec")
for(seed in seed_list) {
    load(paste0('Data/sim_data_', seed, '.rda'))
    true_state = data_format[,'b_true']

    post_prob_S2 = state_results[[seed]][[1]][2, ] / colSums(state_results[[seed]][[1]][1:S,])
    S2_identification = as.numeric(post_prob_S2 >= indiv_select_c)

    # Sensitivity of state 2: Pr(predict S2 | true S2)
    true_S2_ind = S2_identification[true_state == 2]
    sens_S2 = mean(true_S2_ind == 1)

    # Specificity of state 2: Pr(predict not S2 | true not S2)
    true_not_S2_ind = S2_identification[true_state != 2]
    spec_S2 =  mean(true_not_S2_ind == 0)

    sens_and_spec_S2[seed, ] = c(sens_S2, spec_S2)
}

print(paste0("Sensitivity of S2 using c = ", indiv_select_c))
print(summary(sens_and_spec_S2[,"sens"]))

print(paste0("Specificity of S2 using c = ", indiv_select_c))
print(summary(sens_and_spec_S2[,"spec"]))
