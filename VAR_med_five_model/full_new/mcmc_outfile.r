library(latex2exp)
library(tidyverse)
library(gridExtra)

index_seeds = c(1:25)
trialNum = 1
simulation = T

it_num = 2
start_ind = 1

# Parameter initialization -----------------------------------------------------
par_index = list()
par_index$beta = 1:4
par_index$alpha_tilde = 5:20
par_index$upsilon = 21:276
par_index$A = 277:296
par_index$R = 297:312
par_index$zeta = 313:336
par_index$init = 337:340
par_index$omega_tilde = 341:424
par_index$G = 425:440

true_par = rep(0, max(do.call('c', par_index)))
if(simulation) {
    load('Model_out/mcmc_out_1_2it2.rda')
    true_par[par_index$beta] = colMeans(mcmc_out$chain[500:nrow(mcmc_out$chain), mcmc_out$par_index$beta])
    true_par[par_index$alpha_tilde] = colMeans(mcmc_out$chain[500:nrow(mcmc_out$chain), mcmc_out$par_index$alpha_tilde])
    true_par[par_index$upsilon] = colMeans(mcmc_out$chain[500:nrow(mcmc_out$chain), mcmc_out$par_index$upsilon])
    true_par[par_index$A] = colMeans(mcmc_out$chain[500:nrow(mcmc_out$chain), mcmc_out$par_index$A])
    true_par[par_index$R] = colMeans(mcmc_out$chain[500:nrow(mcmc_out$chain), mcmc_out$par_index$R])
    true_par[par_index$zeta] = colMeans(mcmc_out$chain[500:nrow(mcmc_out$chain), mcmc_out$par_index$zeta])
    true_par[par_index$init] = colMeans(mcmc_out$chain[500:nrow(mcmc_out$chain), mcmc_out$par_index$init])
    true_par[par_index$omega_tilde] = colMeans(mcmc_out$chain[500:nrow(mcmc_out$chain), mcmc_out$par_index$omega_tilde])
    
    rm(mcmc_out)    
} else {
    true_par[par_index$beta] = c(0.25, -2, 2, -0.25) # one unit of RBC -> 1 unit increase in hemo in 1 hour
    true_par[par_index$alpha_tilde] = c( -5,   5, -2,  2,
                                         10, -10,  2, -2,
                                         -10,  10,  2, -2,
                                         5,  -5, -2,  2)
    true_par[par_index$upsilon] = c(diag(c(4, 4, 1, 1, 
                                           4, 4, 1, 1, 
                                           4, 4, 1, 1, 
                                           4, 4, 1, 1)))
    true_par[par_index$A] = c(rep(2, 4), rep(-2, 4), rep(0, 4), rep(-2, 4), rep(0, 4))
    true_par[par_index$R] = c(diag(c(9, 9, 9, 9)))
    true_par[par_index$zeta] = c(-3.7405, 2.5, -4.2152,   1, -2.6473,-0.5, -2.1475, -0.2, 
                                 -3.4459,  -1, -2.9404,   1, -3.2151,   1, -3.1778,  1.5, 
                                 -2.0523,   0, -3.4459,-0.2, -3.2404, 2.5, -3.2151,    1)
    true_par[par_index$init] = c(0, 0, 0, 0)
    true_par[par_index$omega_tilde]= 2 * c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                                           -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                                           -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                                           -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
    true_par[par_index$G] = c(diag(c(9, 9, 9, 9)))    
}

load('Data/Dn_omega_names.rda')

labels = c(TeX(r'($\beta_{1}:$ RBC admin effect hemoglobin)'),
           TeX(r'($\beta_{2}:$ RBC admin effect heart rate)'),
           TeX(r'($\beta_{3}:$ RBC admin effect MAP)'),
           TeX(r'($\beta_{4}:$ RBC admin effect lactate)'),
           TeX(r'($\tilde{\alpha}_{1,1}:$ state 2 slope hemoglobin)'),
           TeX(r'($\tilde{\alpha}_{2,1}:$ state 3 slope hemoglobin)'),
           TeX(r'($\tilde{\alpha}_{3,1}:$ state 4 slope hemoglobin)'),
           TeX(r'($\tilde{\alpha}_{4,1}:$ state 5 slope hemoglobin)'),
           TeX(r'($\tilde{\alpha}_{1,2}:$ state 2 slope heart rate)'),
           TeX(r'($\tilde{\alpha}_{2,2}:$ state 3 slope heart rate)'),
           TeX(r'($\tilde{\alpha}_{3,2}:$ state 4 slope heart rate)'),
           TeX(r'($\tilde{\alpha}_{4,2}:$ state 5 slope heart rate)'),
           TeX(r'($\tilde{\alpha}_{1,3}:$ state 2 slope MAP)'),
           TeX(r'($\tilde{\alpha}_{2,3}:$ state 3 slope MAP)'),
           TeX(r'($\tilde{\alpha}_{3,3}:$ state 4 slope MAP)'),
           TeX(r'($\tilde{\alpha}_{4,3}:$ state 5 slope MAP)'),
           TeX(r'($\tilde{\alpha}_{1,4}:$ state 2 slope lactate)'),
           TeX(r'($\tilde{\alpha}_{2,4}:$ state 3 slope lactate)'),
           TeX(r'($\tilde{\alpha}_{3,4}:$ state 4 slope lactate)'),
           TeX(r'($\tilde{\alpha}_{4,4}:$ state 5 slope lactate)'),
           paste0(TeX(r'($\Upsilon_{\alpha}$( )'), 1:16, ", ", rep(1:16, each = 16), " )"), 
           TeX(r'(logit $A_{1}(hemoglobin)$)'), TeX(r'(logit $A_{1}(heart rate)$)'),
           TeX(r'(logit $A_{1}(MAP)$)'), TeX(r'(logit $A_{1}(lactate)$)'),
           TeX(r'(logit $A_{2}(hemoglobin)$)'), TeX(r'(logit $A_{2}(heart rate)$)'),
           TeX(r'(logit $A_{2}(MAP)$)'), TeX(r'(logit $A_{2}(lactate)$)'),
           TeX(r'(logit $A_{3}(hemoglobin)$)'), TeX(r'(logit $A_{3}(heart rate)$)'),
           TeX(r'(logit $A_{3}(MAP)$)'), TeX(r'(logit $A_{3}(lactate)$)'),
           TeX(r'(logit $A_{4}(hemoglobin)$)'), TeX(r'(logit $A_{4}(heart rate)$)'),
           TeX(r'(logit $A_{4}(MAP)$)'), TeX(r'(logit $A_{4}(lactate)$)'),
           TeX(r'(logit $A_{5}(hemoglobin)$)'), TeX(r'(logit $A_{5}(heart rate)$)'),
           TeX(r'(logit $A_{5}(MAP)$)'), TeX(r'(logit $A_{5}(lactate)$)'),
           TeX(r'($R_{1,1}:$ Var(hemo) )'), TeX(r'($R_{2,1}:$ Cov(hemo, hr) )'),
           TeX(r'($R_{3,1}:$ Cov(hemo, map) )'), TeX(r'($R_{4,1}:$ Cov(hemo, lact) )'),
           TeX(r'($R_{1,2}:$ Cov(hr, hemo) )'), TeX(r'($R_{2,2}:$ Var(hr) )'),
           TeX(r'($R_{3,2}:$ Cov(hr, map) )'), TeX(r'($R_{4,2}:$ Cov(hr, lact) )'),
           TeX(r'($R_{1,3}:$ Cov(map, hemo) )'), TeX(r'($R_{2,3}:$ Cov(map, hr) )'),
           TeX(r'($R_{3,3}:$ Var(map) )'), TeX(r'($R_{4,3}:$ Cov(map, lact) )'),
           TeX(r'($R_{1,3}:$ Cov(lact, hemo) )'), TeX(r'($R_{2,3}:$ Cov(lact, hr) )'),
           TeX(r'($R_{3,3}:$ Cov(lact, map) )'), TeX(r'($R_{4,3}:$ Var(lact) )'),
           TeX(r'($\zeta_{1,1}:$ S1 --> S2)'), TeX(r'($\zeta_{2,1}:$ S1 --> S2)'),
           TeX(r'($\zeta_{1,2}:$ S1 --> S4)'), TeX(r'($\zeta_{2,2}:$ S1 --> S4)'),
           TeX(r'($\zeta_{1,3}:$ S2 --> S3)'), TeX(r'($\zeta_{2,3}:$ S2 --> S3)'),
           TeX(r'($\zeta_{1,4}:$ S2 --> S4)'), TeX(r'($\zeta_{2,4}:$ S2 --> S4)'),
           TeX(r'($\zeta_{1,5}:$ S3 --> S1)'), TeX(r'($\zeta_{2,5}:$ S3 --> S1)'),
           TeX(r'($\zeta_{1,6}:$ S3 --> S2)'), TeX(r'($\zeta_{2,6}:$ S3 --> S2)'),
           TeX(r'($\zeta_{1,7}:$ S3 --> S4)'), TeX(r'($\zeta_{2,7}:$ S3 --> S4)'),
           TeX(r'($\zeta_{1,8}:$ S4 --> S2)'), TeX(r'($\zeta_{2,8}:$ S4 --> S2)'),
           TeX(r'($\zeta_{1,9}:$ S4 --> S5)'), TeX(r'($\zeta_{2,9}:$ S4 --> S5)'),
           TeX(r'($\zeta_{1,10}:$ S5 --> S1)'), TeX(r'($\zeta_{2,10}:$ S5 --> S1)'),
           TeX(r'($\zeta_{1,11}:$ S5 --> S2)'), TeX(r'($\zeta_{2,11}:$ S5 --> S2)'),
           TeX(r'($\zeta_{1,12}:$ S5 --> S4)'), TeX(r'($\zeta_{2,12}:$ S5 --> S4)'),
           "logit Pr(init S2)", "logit Pr(init S3)",
           "logit Pr(init S4)", "logit Pr(init S5)",
           paste0("mean (hr): ", Dn_omega_names[1:34]), 
           paste0("mean (map): ", Dn_omega_names[35:84]), 
           "G(hemo)", "G(hemo, hr)", "G(hemo, map)", "G(hemo, lact)", 
           "G(hr, hemo)", "G(hr)", "G(hr, map)", "G(hr, lact)",
           "G(map, hemo)", "G(map, hr)", "G(map)", "G(map, lact)",
           "G(lact, hemo)", "G(lact, hr)", "G(lact, map)", "G(lact)")

direction_med = c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                  -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                  -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                  -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)

if(simulation) {
    # Estimate the initial state probabilities
    init_prob_mat = matrix(nrow = length(index_seeds), ncol = 5)
    c_s = 1
    for(seed in index_seeds) {
        load(paste0('Data/sim_data_', seed, '.rda'))
        first_ind = c(0, which(diff(data_format[,"EID"]) != 0)) + 1
        init_state = data_format[first_ind, "b_true"]
        init_prob_mat[c_s, ] = c(sum(init_state == 1), sum(init_state == 2),
                                 sum(init_state == 3), sum(init_state == 4), 
                                 sum(init_state == 5)) / length(init_state)
        c_s = c_s + 1
    }
    
    init_par_est = colMeans(init_prob_mat)
    true_par[par_index$init[1]] = log(init_par_est[2] / (1 - sum(init_par_est[2:5])))
    true_par[par_index$init[2]] = log(init_par_est[3] / (1 - sum(init_par_est[2:5])))
    true_par[par_index$init[3]] = log(init_par_est[4] / (1 - sum(init_par_est[2:5])))
    true_par[par_index$init[4]] = log(init_par_est[5] / (1 - sum(init_par_est[2:5])))
}

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------

chain_list = list()
A_list = list()
par_means = vector(mode = 'list', length = length(index_seeds))
ind = 0
post_median_mat = NULL
covg = NULL

for(seed in index_seeds){
    
    it_seq = 1:(it_num - start_ind + 1)
    covg_val = FALSE
    
    if(simulation) {
        check_name = paste0('Model_out/mcmc_out_',trialNum, '_', seed, 
                            'it', it_num,'_sim.rda')    
    } else {
        check_name = paste0('Model_out/mcmc_out_',trialNum, '_', seed, 
                            'it', it_num,'.rda')
    }
    
    if(file.exists(check_name)) {
        
        for(it in it_seq) {
            
            if(simulation) {
                load(paste0('Model_out/mcmc_out_',trialNum, '_', seed, 
                            'it', it + (start_ind - 1),'_sim.rda'))    
            } else {
                load(paste0('Model_out/mcmc_out_',trialNum, '_', seed, 
                            'it', it + (start_ind - 1),'.rda'))
            }
            
            print(paste0(seed, ": ", it))
            print("accept")
            print(mcmc_out$accept)
            
            par_index = mcmc_out$par_index
            
            if(it == 1) {
                ind = ind + 1
                chain_list[[ind]] = mcmc_out$chain[500:1000, ]
                A_list[[ind]] = t(do.call(cbind, mcmc_out$alpha_i))
            } else {
                chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out$chain)
                A_list[[ind]] = rbind(A_list[[ind]], t(do.call(cbind, mcmc_out$alpha_i)))
            }
            
            rm(mcmc_out)
        }
        
        covg_val = TRUE
        
    } else {
        print(paste0("Missing! ", check_name))
    }
    
    if(covg_val) {
        post_median_temp = matrix(apply(chain_list[[ind]], 2, median), nrow = 1)
        post_median_mat = rbind(post_median_mat, post_median_temp)
        
        par_means[[seed]] = colMeans(chain_list[[ind]])
        
        covg_low = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.025)})
        covg_high = apply(chain_list[[ind]], 2, function(x){quantile(x, prob=.975)})
        covg_temp = as.numeric(true_par <= covg_high & true_par >= covg_low)
        covg = rbind(covg, covg_temp)    
    }
}

save(par_means, file = paste0('Model_out/par_means_it', it_num, '_', as.numeric(simulation), '.rda'))

upsilon_ind = matrix(1:length(par_index$upsilon), ncol = 16)
upsilon_ind[upper.tri(upsilon_ind, diag = F)] = 0
upsilon_ind = c(upsilon_ind)
upsilon_ind = upsilon_ind[upsilon_ind != 0]

stacked_chains = do.call( rbind, chain_list)

post_med_stack = c(apply(stacked_chains, 2, median))
save(post_med_stack, file = paste0('Model_out/post_med_stack_it', it_num, '_', as.numeric(simulation), '.rda'))

pdf(paste0('Plots/trace_plot_', as.numeric(simulation), '_it', it_num, '.pdf'))
if(simulation) {
    par(mfrow=c(4, 2))    
} else {
    par(mfrow=c(3, 2))
}

lab_ind = 0
med_print_save = NULL
for(s in names(par_index)){
    temp_par = par_index[[s]]
    if(s == "upsilon") {
        temp_par = temp_par[upsilon_ind]
    }
    
    for(r in temp_par){
        # lab_ind = lab_ind + 1
        lab_ind = r
        parMean = round( mean(stacked_chains[,r], na.rm = T), 4)
        parMedian = round( median(stacked_chains[,r], na.rm = T), 4)
        upper = quantile( stacked_chains[,r], prob=.975, na.rm = T)
        lower = quantile( stacked_chains[,r], prob=.025, na.rm = T)
        
        if(s == "A") {
            med_trans = exp(parMedian) / (1 + exp(parMedian))
            low_trans = exp(lower) / (1 + exp(lower))
            upp_trans = exp(upper) / (1 + exp(upper))
            print(paste0(labels[lab_ind], ": ", round(med_trans, 3), " [", round(low_trans, 3), ", ", round(upp_trans, 3), "]"))
        } else {
            if(s == "omega_tilde") {
                if(0 < lower || 0 > upper) {
                    if(parMedian * direction_med[r - 340] < 0) {
                        omega_label = paste0("\\textcolor{red}{$", round(parMedian, 3), "^*$}")
                    } else {
                        omega_label = paste0("$", round(parMedian, 3), "^*$")   
                    }
                } else {
                    omega_label = paste0("$", round(parMedian, 3), "$")
                }
                
                med_print_save = c(med_print_save, omega_label)
            }
            print(paste0(labels[lab_ind], ": ", round(parMedian, 3), " [", round(lower, 3), ", ", round(upper, 3), "]"))    
        }
        
        y_limit = range(stacked_chains[,r])
        plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,nrow(chain_list[[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"))
        
        for(seed in 1:length(chain_list)) { lines( chain_list[[seed]][,r], type='l', col=seed) }
        
        if(simulation) {
            x_label = paste0('Mean =',toString(parMean),
                             ' Median =',toString(parMedian),
                             ' True =', round(true_par[r], 3)) 
        } else {
            x_label = paste0('Mean =',toString(parMean),
                             ' Median =',toString(parMedian))     
        }
        
        hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, freq=FALSE,
              xlab=x_label, main = "")
        abline( v=upper, col='red', lwd=2, lty=2)
        abline( v=lower, col='purple', lwd=2, lty=2)
        if(simulation) { abline( v=true_par[r], col='green', lwd=2, lty=2) }
    }   
}
dev.off()

# Plot the sampled alpha_i
alpha_i_names = c(TeX(r'($[\alpha^{(i)}_{*}]_{.,1}:$ hemoglobin slopes)'),
                  TeX(r'($[\alpha^{(i)}_{*}]_{.,2}:$ heart rate slopes)'),
                  TeX(r'($[\alpha^{(i)}_{*}]_{.,3}:$ MAP slopes)'),
                  TeX(r'($[\alpha^{(i)}_{*}]_{.,4}:$ lactate slopes)'))

all_a = do.call(rbind, A_list)
VP <- vector(mode="list", length = length(alpha_i_names))

for(a in 1:4) {
    state_per_a = (4*(a-1) + 1):(4*a)
    
    yVar = disc_type = NULL
    
    yVar = c(all_a[,state_per_a])
    
    print(summary(all_a[,state_per_a[1]]))
    
    disc_type = c(rep("state 2", nrow(all_a)), rep("state 3", nrow(all_a)),
                  rep("state 4", nrow(all_a)), rep("state 5", nrow(all_a)))
    
    plot_df = data.frame(yVar = yVar, disc_type = disc_type)
    plot_df$disc_type <- factor(plot_df$disc_type, levels = unique(plot_df$disc_type))
    VP[[a]] = ggplot(plot_df, aes(x=disc_type, y = yVar, fill = disc_type)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1) +
        ggtitle(alpha_i_names[a]) +
        ylab("") + #paste0("Parameter Value: ", round(true_par[r], 3))
        xlab(TeX(r'(sampled $\alpha^{(i)}_{*}$)')) +
        theme(legend.position = "none", 
              panel.background = element_rect(fill = "white", colour = "grey", 
                                              linetype = 'solid', linewidth = 1),
              panel.grid.major = element_line(linetype = 'solid',
                                              colour = "grey") )
    
}

pdf(paste0('Plots/sampled_alpha_', as.numeric(simulation), '_it', it_num, '.pdf'))
grid.arrange(VP[[1]], VP[[2]], nrow = 2)
grid.arrange(VP[[3]], VP[[4]], nrow = 2)
dev.off()


# Print Medication information in a clear format
direction_text = rep(" ", length(direction_med))
direction_text[direction_med==1] = "$\\uparrow$"
direction_text[direction_med==-1] = "$\\downarrow$"
med_key = read.csv('Data/Med_chart.csv', na.strings = "")
colnames(med_key) = c('id', 'hr', 'map', 'onset', 'offset', 'time_check', 'X1')
med_key$id = paste0(" ", med_key$id)
med_key$id[med_key$id == " Metoprolol IV"] = " METOPROLOL"
library(stringr)

all_med_names = med_key$id

all_med_names <- str_trim(all_med_names, "right")

count_med_effects = 0
for(mn in 1:length(all_med_names)) {
    check_med = grepl(all_med_names[mn], Dn_omega_names, fixed = TRUE)
    if(sum(check_med) > 0) {
    
        dn_loc = which(check_med)    
        dn_name = Dn_omega_names[dn_loc]
        
        # order: hr cont, hr disc, map cont, map disc
        info_print = c("--", "--", "--", "--")
        dir_print  = c("", "", "", "")
        for(l in 1:length(dn_loc)) {
            if(dn_loc[l] < 35) {
                if(grepl('1', dn_name[l], fixed=TRUE)) {
                    info_print[1] = med_print_save[dn_loc[l]]
                    dir_print[1] = direction_text[dn_loc[l]]
                } else if(grepl('0', dn_name[l], fixed=TRUE)) {
                    info_print[2] = med_print_save[dn_loc[l]]
                    dir_print[2] = direction_text[dn_loc[l]]
                }
            } else {
                if(grepl('1', dn_name[l], fixed=TRUE)) {
                    info_print[3] = med_print_save[dn_loc[l]]
                    dir_print[3] = direction_text[dn_loc[l]]
                } else if(grepl('0', dn_name[l], fixed=TRUE)) {
                    info_print[4] = med_print_save[dn_loc[l]]
                    dir_print[4] = direction_text[dn_loc[l]]
                }
            }
        }
        
        count_med_effects = count_med_effects + sum(info_print != "--")
        
        cat(all_med_names[mn], " & ", dir_print[1], dir_print[2], " & ", info_print[1], " & ", info_print[2],
                     " & ", dir_print[3], dir_print[4], " & ", info_print[3], " & ", info_print[4]," \\\\", '\n')
        
    } 
}


print(paste0("Total med effects: ", count_med_effects))

if(simulation) {
    pdf(paste0('Plots/box_plots_it', it_num, '_sim.pdf'))
    par(mfrow=c(4, 4), mar = c(2, 2, 3, 1))
    lab_ind = 0
    for(s in names(par_index)){
        temp_par = par_index[[s]]
        if(s == "upsilon") {
            temp_par = temp_par[upsilon_ind]
        }
        for(r in temp_par){
            
            lab_ind = r
            boxplot(post_median_mat[,r], main = labels[r],
                    xlab = " ",
                    ylab = " ",
                    xaxs = "i",
                    yaxs = "i")
            title(xlab=paste0('95% Covg = ', round(mean(covg[,r]), 4), 
                           '\n', 'True value = ', round(true_par[r], 4)), 
                  line=1)
            abline(h = true_par[r], col = 'red')
        }
    }
    dev.off()
    
    pdf(paste0('Plots/box_plots_alpha_it', it_num, '_sim.pdf'))
    par(mfrow=c(4, 4), mar = c(2, 2, 3, 1))
    
    temp_par = par_index[[2]]
    
    for(r in temp_par){
        
        lab_ind = r
        boxplot(post_median_mat[,r], main = labels[r],
                xlab = " ",
                ylab = " ",
                xaxs = "i",
                yaxs = "i")
        title(xlab=paste0('95% Covg = ', round(mean(covg[,r]), 4), 
                          '\n', 'True value = ', round(true_par[r], 4)), line=1)
        abline(h = true_par[r], col = 'red')
    }
    dev.off()
}




