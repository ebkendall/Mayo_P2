library(matrixStats)
library(plotrix)
library(splines)

simulation = T

trialNum = 1
args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])

states_per_step = 0
steps_per_it = 1
S = 5

it_num = 1
it_seq = 1:it_num

seed_list = c(1:50)

# Mode of the state sequences -------------------------------------------------
Mode <- function(x) {
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
}

# Load the model output -------------------------------------------------------
state_results = list()

ind = 0
flag_ind = -1
for(seed in seed_list){
    
    B_chain   = NULL
    
    for(it in it_seq) {
        
        file_name = paste0('Model_out/mcmc_out_', seed, '_', seed, 'it', 
                           it, '_samp', sampling_num, '_', states_per_step,
                           '_', steps_per_it,'_sim.rda') 
        
        if(file.exists(file_name)) {
            load(file_name)
            print(paste0(ind+1, ": ", file_name))
            # keep_ind = seq(1, nrow(mcmc_out_temp$B_chain), by = 5)
            
            # if(it == 1) {
                B_chain   = mcmc_out_temp$B_chain[1000:2000, ]
            # } else {
            #     B_chain   = rbind(B_chain, mcmc_out_temp$B_chain[keep_ind, ])
            # }
            
            rm(mcmc_out_temp)    
        } else {
            print(paste0("Missing: ", file_name))
        }
    }
    
    ind = ind + 1
    
    if(!is.null(B_chain)) {
        state_results[[ind]] = matrix(nrow = S+1, ncol = ncol(B_chain))
        for(jj in 1:S) {
            state_results[[ind]][jj, ] = apply(B_chain, 2, function(x,jj){sum(x == jj)}, jj)
        }
        state_results[[ind]][S+1, ] = apply(B_chain, 2, Mode)
    } else {
        state_results[[ind]] = NA
        flag_ind = c(flag_ind, ind)
    }
}

# Computing ROC curve information ---------------------------------------------
mode_correct = rep(NA, length(seed_list))
c = seq(0, 1, by = 0.001)
chosen_c = NULL

pdf(paste0('Plots/ROC_', trialNum, '_samp', sampling_num, '_',
           states_per_step, '_', steps_per_it, '.pdf'))
par(mfrow=c(3, 2))
for(seed in seed_list) {
    print(seed)
    if(!(seed %in% flag_ind)) {
        load(paste0('Data_sim/use_data_', seed, '.rda'))
        true_state = use_data[,'b_true']
        
        # Quick summary 
        mode_correct[seed] = mean(state_results[[seed]][S+1, ] == true_state)
        
        # State 2 focus
        state_2_binary = as.numeric(true_state == 2)
        state_1345_binary = as.numeric(true_state != 2)
        
        cumulative_post_prob = state_results[[seed]][2, ] / colSums(state_results[[seed]][1:S,])
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
        
        slope_mat = cbind(diff(c_results$true_pos), diff(c_results$false_pos))
        slope_mat = cbind(c[1:nrow(slope_mat)], slope_mat)
        colnames(slope_mat) = c("c", "delta_y", "delta_x")
        
        slope_mat = slope_mat[slope_mat[,"delta_x"] != 0, ]
        slopes = slope_mat[,"delta_y"] / slope_mat[,"delta_x"]
        slope_mat = cbind(slope_mat, slopes)
        
        approx_1_slope = slope_mat[order(abs(slope_mat[,"slopes"] - 1))[1], "slopes"]
        c_slope_1 = slope_mat[slope_mat[,"slopes"] == approx_1_slope, "c"]
        
        chosen_c = rbind(chosen_c, c_results[c_results$c %in% c_slope_1, ])
        
        # AUC (estimated using Trapezoid rule)
        area = 0
        for(k in 2:nrow(c_results)) {
            area = area + 0.5 * (c_results$true_pos[k] + c_results$true_pos[k-1]) * 
                (c_results$false_pos[k] - c_results$false_pos[k-1])
        }
        
        plot(c_results$false_pos, c_results$true_pos, xlim = c(0,1),
             xlab = "FPR", ylab = "TPR", 
             main = paste0(seed, ", State 2 AUC = ", round(area, digits = 4)))
        points(c_results$false_pos[c_results$c %in% c_slope_1],
               c_results$true_pos[c_results$c %in% c_slope_1],
               col = 'red', cex = 2)
        
    }
}
dev.off()

print("Summary of state identification based on posterior mode: ")
print(summary(mode_correct))

print("Summary of the calibrated minimum probability")
print(summary(chosen_c$c))

indiv_select_c = mean(chosen_c$c)
sens_and_spec_S2 = matrix(nrow = length(seed_list), ncol = 2)
colnames(sens_and_spec_S2) = c("sens", "spec")
for(seed in seed_list) {
    if(!(seed %in% flag_ind)) {
        load(paste0('Data_sim/use_data_', seed, '.rda'))
        true_state = use_data[,'b_true']
        
        post_prob_S2 = state_results[[seed]][2, ] / colSums(state_results[[seed]][1:S,])
        S2_identification = as.numeric(post_prob_S2 >= indiv_select_c)
        
        # Sensitivity of state 2: Pr(predict S2 | true S2)
        true_S2_ind = S2_identification[true_state == 2]
        sens_S2 = mean(true_S2_ind == 1)
        
        # Specificity of state 2: Pr(predict not S2 | true not S2)
        true_not_S2_ind = S2_identification[true_state != 2]
        spec_S2 =  mean(true_not_S2_ind == 0)
        
        sens_and_spec_S2[seed, ] = c(sens_S2, spec_S2)
    }
}

print(paste0("Sensitivity of S2 using c = ", indiv_select_c))
print(summary(sens_and_spec_S2[,"sens"]))

print(paste0("Specificity of S2 using c = ", indiv_select_c))
print(summary(sens_and_spec_S2[,"spec"]))
