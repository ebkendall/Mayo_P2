args = commandArgs(TRUE)
before_t1 = as.numeric(args[1])

index_seeds = c(1:25)
it_num = 1
S = 3

# Mode of the state sequences -------------------------------------------------
Mode <- function(x) {
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
}

# Load the model output -------------------------------------------------------
state_results = vector(mode = 'list', length = length(index_seeds))

ind = 0
for(seed in index_seeds){
    
    it_seq = 1:it_num
    
    for(it in it_seq) {
        file_name = paste0('Model_out/mcmc_out_',seed, '_it_', it, '_', before_t1,'.rda')
        if(file.exists(file_name)) {
            load(file_name)
            print(paste0(seed, ": ", file_name))
            
            par_index = mcmc_out$par_index
            
            if(it == 1) {
                ind = ind + 1
                B_chain = mcmc_out$B_chain[200:1000, ]
            } else {
                B_chain = rbind(B_chain, mcmc_out_temp$B_chain)
            }
            
            rm(mcmc_out)
        } else {
            print(paste0("Missing! ", file_name))
        }
    }
    
    state_results[[seed]] = matrix(nrow = S+1, ncol = ncol(B_chain))
    for(jj in 1:S) {
        state_results[[seed]][jj, ] = apply(B_chain, 2, function(x,jj){sum(x == jj)}, jj)
    }
    state_results[[seed]][S+1, ] = apply(B_chain, 2, Mode)
    
    rm(B_chain)
}

# Summary stats of state identification ----------------------------------------
for(seed in 1:length(index_seeds)) {
    
    load(paste0('Data/data_format', seed, '.rda'))
    ss_true = data_format[,"state"]
    
    if(length(ss_true) != length(state_results[[seed]][S+1, ])) {
        print(paste0("ERROR (seed ", seed, ")"))
    } else {
        prop_all = sum(ss_true == state_results[[seed]][S+1, ]) / length(state_results[[seed]][S+1, ])
        print(paste0("Across all subjects (seed ", seed, ") = ", prop_all))
    }
    
    EIDs = unique(data_format[,"id"])
    
    prop_sub = rep(0, length(EIDs))
    
    for(j in 1:length(EIDs)) {
        sub_ind_j = which(data_format[,"id"] == EIDs[j])
        ss_true_j = ss_true[sub_ind_j]
        state_seq_mode_j = state_results[[seed]][S+1, sub_ind_j]
        
        prop_sub[j] = sum(ss_true_j == state_seq_mode_j) / length(state_seq_mode_j)
    }
    
    print(paste0("Subject specific (seed ", seed, ")"))
    print(summary(prop_sub))
    
    # Sensitivity of state 2: Pr(predict S2 | true S2)
    predict_at_true_S2 = state_results[[seed]][S+1, (ss_true == 2)]
    print(paste0("Sensitivity of S2 = ", mean(predict_at_true_S2 == 2)))
    
    # Specificity of state 2: Pr(predict not S2 | true not S2)
    predict_not_S2 =  state_results[[seed]][S+1, (ss_true != 2)]
    print(paste0("Specificity of S2 = ", mean(predict_not_S2 != 2)))
}