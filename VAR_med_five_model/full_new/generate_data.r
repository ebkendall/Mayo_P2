library(mvtnorm, quietly=T)

# Parameter initialization -----------------------------------------------------
par_index = list()
par_index$beta = 1:4
par_index$alpha_tilde = 5:24
par_index$upsilon = 25:424
par_index$A = 425:444
par_index$R = 445:460
par_index$zeta = 461:484
par_index$init = 485:488
par_index$omega_tilde = 489:572

# Initializing true parameter values -------------------------------------------
par = rep(0, max(do.call('c', par_index)))

par[par_index$beta] = c(0.25, -2, 2, -0.25) # one unit of RBC -> 1 unit increase in hemo in 1 hour
par[par_index$alpha_tilde] = c( 50,  -5,   5, -2,  2,
                               100,  10, -10,  2, -2,
                               100, -10,  10,  2, -2,
                                50,   5,  -5, -2,  2)
par[par_index$upsilon] = c(diag(c(25, 4, 4, 1, 1, 
                                  25, 4, 4, 1, 1, 
                                  25, 4, 4, 1, 1, 
                                  25, 4, 4, 1, 1)))
par[par_index$A] = c(rep(2, 4), rep(-2, 4), rep(0, 4), rep(-2, 4), rep(0, 4))
par[par_index$R] = c(diag(c(9, 9, 9, 9)))
#    transitions:          1->2,         1->4,         2->3,         2->4, 
#                          3->1,         3->2,         3->4,         4->2, 
#                          4->5,         5->1,         5->2,         5->4
par[par_index$zeta] = c(-3.7405, 2.5, -4.2152,   1, -2.6473,-0.5, -2.1475, -0.2, 
                        -3.4459,  -1, -2.9404,   1, -3.2151,   1, -3.1778,  1.5, 
                        -2.0523,   0, -3.4459,-0.2, -3.2404, 2.5, -3.2151,    1)
par[par_index$init] = c(0, 0, 0, 0)
par[par_index$omega_tilde]= 2 * c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                                  -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                                  -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                                  -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)

n_state = 5

# Parameter initialization -----------------------------------------------------
beta = par[par_index$beta]
alpha_tilde = matrix(par[par_index$alpha_tilde], ncol = 4)
Upsilon = matrix(par[par_index$upsilon], ncol = 20)
vec_A = par[par_index$A]
A_scaled = exp(vec_A) / (1 + exp(vec_A)) # Support (0,1)
A_total = matrix(A_scaled, ncol = 5)

# columns: hemo, hr, map, lactate
R = matrix(par[par_index$R], ncol = 4)

zeta = matrix(par[par_index$zeta], nrow = 2)
colnames(zeta) = c('(1) 1->2', '(2) 1->4','(3) 2->3', '(4) 2->4', '(5) 3->1', 
                   '(6) 3->2', '(7) 3->4','(8) 4->2', '(9) 4->5', '(10) 5->1', 
                   '(11) 5->2', '(12) 5->4')

init_logit = par[par_index$init]
init_logit = c(0, init_logit)

omega_tilde = par[par_index$omega_tilde]
# ------------------------------------------------------------------------------

# for(df_num in 1:50) {
args = commandArgs(TRUE)
df_num = as.numeric(args[1])

print(paste0("Df ", df_num))
set.seed(df_num)

load('Data/data_format_train_update.rda')
load('Data/Dn_omega_update.rda')
load('Data/Dn_omega_names.rda')

EIDs = unique(data_format[,"EID"])
N = length(unique(data_format[,"EID"]))
Y = data_format[, c('EID', 'hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
x = data_format[,c('n_RBC_admin'), drop=F]
z = cbind(1, data_format[,c('RBC_ordered'), drop=F])

otype = !is.na(Y[, c('hemo','hr','map','lactate')])
colnames(otype) = c('hemo','hr','map','lactate')

# Clinic Rule information --------------------------------------------------
non_zero_clinic = unique(data_format[data_format[,"clinic_rule"] != 0, "EID"])
clinic_rbc_combos = matrix(0, ncol = 3, nrow = length(non_zero_clinic))
clinic_rbc_combos[,1] = non_zero_clinic
for(c in 1:length(non_zero_clinic)) {
    clinic_rbc_combos[c,2:3] = c(unique(data_format[data_format[,"EID"] == non_zero_clinic[c], "RBC_rule"]),
                                 unique(data_format[data_format[,"EID"] == non_zero_clinic[c], "clinic_rule"]))
}

# Indicator variable of the first RBC to indicate bleed event --------------
bleed_pat = unique(data_format[data_format[,"RBC_rule"] != 0, "EID"])
bleed_indicator = rep(0, nrow(data_format))
for(i in 1:length(bleed_pat)) {
    sub_dat = data_format[data_format[,"EID"] == bleed_pat[i], ]
    
    # Check in any 12 hour period
    max_time = tail(sub_dat[,"time"], 1)
    when_rbc = c(1, which(diff(sub_dat[,"n_RBC_admin"]) != 0))
    
    for(j in 1:length(when_rbc)) {
        s_time = sub_dat[when_rbc[j], "time"]
        e_time_12 = s_time + 720
        RBC_diff_12 = 0
        
        if (e_time_12 <= max_time) {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            ind_12 = order(abs(sub_dat[,"time"] - e_time_12))[1]
            RBC_diff_12 = sub_dat[ind_12, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        } else {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
            RBC_diff_12 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        }
        
        if(RBC_diff_12 >=3) {
            
            admin_times = sub_dat[sub_dat[,"RBC_admin"] != 0, "time"]
            
            a_t = which(admin_times >= s_time & admin_times < e_time_12)
            
            # first_time = admin_times[a_t[1]]       # before first of three RBCs
            first_time = admin_times[tail(a_t, 1)] # before last of three RBCs
            
            order_times = sub_dat[sub_dat[,"RBC_ordered"] != 0, "time"]
            if(sum(order_times <= first_time) == 0) {
                first_order_time = first_time
            } else {
                first_order_time = max(order_times[order_times <= first_time])   
            }
            
            bleed_indicator[data_format[,"EID"] == bleed_pat[i] & 
                                data_format[,"time"] == first_order_time] = 1
            break
        }
    }
}

rm(data_format)
# --------------------------------------------------------------------------

alpha_i_mat = vector(mode = "list", length = N)
Dn_omega_sim = vector(mode = 'list', length = N)

data_format = NULL
rbc_bleed_correct = NULL
initial_state_vec = NULL

for(i in 1:N) {
    
    id_num = EIDs[i]
    
    Dn_omega_sim[[i]] = Dn_omega[[i]]
    D_i_omega = Dn_omega_sim[[i]]
    
    rbc_rule = as.logical(head(Y[Y[,'EID']==as.numeric(id_num),"RBC_rule"], 1))
    correct_bleed = T
    if(rbc_rule) {correct_bleed = F}
    
    n_i = sum(Y[,'EID']==as.numeric(id_num))
    
    # m_i = n_i + rpois(n = 1, lambda = 50)
    m_i = n_i
    
    x_i = x[ Y[,'EID']==as.numeric(id_num),, drop=F]
    z_i = z[ Y[,'EID']==as.numeric(id_num),, drop=F]
    
    bleed_ind_i = bleed_indicator[Y[,'EID']==as.numeric(id_num)]
    
    if(length(D_i_omega) != n_i) {print(paste0("issue n_i: ", i))}
    
    # Sample a long string of states ---------------------------------------
    big_b_i = rep(NA, m_i)
    big_b_i[1] = 1
    for(k in 2:m_i) {
        if(k <= m_i - n_i) {
            q1   = exp(zeta[1, 1]) 
            q2   = exp(zeta[1, 2])
            q3   = exp(zeta[1, 3])
            q4   = exp(zeta[1, 4])
            q5   = exp(zeta[1, 5]) 
            q6   = exp(zeta[1, 6])
            q7   = exp(zeta[1, 7])
            q8   = exp(zeta[1, 8])
            q9   = exp(zeta[1, 9]) 
            q10  = exp(zeta[1, 10])
            q11  = exp(zeta[1, 11])
            q12  = exp(zeta[1, 12])
            
            Q = matrix(c(  1,   q1,  0,  q2,  0,
                           0,    1, q3,  q4,  0,
                           q5,   q6,  1,  q7,  0,
                           0,   q8,  0,   1, q9,
                           q10,  q11,  0, q12,  1), ncol=5, byrow=T)
            
            P_i = Q / rowSums(Q)
            
            big_b_i[k] = sample(1:5, size=1, prob=P_i[big_b_i[k-1],])
        } else {
            q1   = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  1, drop=F]) 
            q2   = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  2, drop=F])
            q3   = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  3, drop=F])
            q4   = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  4, drop=F])
            q5   = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  5, drop=F]) 
            q6   = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  6, drop=F])
            q7   = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  7, drop=F])
            q8   = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  8, drop=F])
            q9   = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  9, drop=F]) 
            q10  = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  10, drop=F])
            q11  = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  11, drop=F])
            q12  = exp(z_i[k - (m_i - n_i),, drop=F] %*% zeta[,  12, drop=F])
            
            Q = matrix(c(   1,   q1,  0,  q2,  0,
                            0,    1, q3,  q4,  0,
                            q5,   q6,  1,  q7,  0,
                            0,   q8,  0,   1, q9,
                            q10,  q11,  0, q12,  1), ncol=5, byrow=T)
            
            P_i = Q / rowSums(Q)
            
            big_b_i[k] = sample(1:5, size=1, prob=P_i[big_b_i[k-1],])
        }
    }
    
    # Select n_i of the states for the true latent state sequence ----------
    b_i = tail(big_b_i, n_i)
    before_b_i = big_b_i[1:(m_i - n_i + 1)]
    
    # Generate data --------------------------------------------------------
    Y_i = matrix(nrow = n_i, ncol = 4)
    D_i = vector(mode = 'list', length = n_i)
    X_i = vector(mode = 'list', length = n_i)
    
    vec_alpha_i = rmvnorm( n=1, mean=c(alpha_tilde), sigma=Upsilon)
    
    alpha_i_mat[[i]] = matrix(vec_alpha_i, ncol = 1)
    
    A_init = A_total[,b_i[1]]
    Gamma = matrix(c(R[1,1] / (1 - A_init[1] * A_init[1]), 
                     R[1,2] / (1 - A_init[1] * A_init[2]),
                     R[1,3] / (1 - A_init[1] * A_init[3]), 
                     R[1,4] / (1 - A_init[1] * A_init[4]),
                     R[2,1] / (1 - A_init[2] * A_init[1]), 
                     R[2,2] / (1 - A_init[2] * A_init[2]),
                     R[2,3] / (1 - A_init[2] * A_init[3]), 
                     R[2,4] / (1 - A_init[2] * A_init[4]),
                     R[3,1] / (1 - A_init[3] * A_init[1]), 
                     R[3,2] / (1 - A_init[3] * A_init[2]),
                     R[3,3] / (1 - A_init[3] * A_init[3]), 
                     R[3,4] / (1 - A_init[3] * A_init[4]),
                     R[4,1] / (1 - A_init[4] * A_init[1]), 
                     R[4,2] / (1 - A_init[4] * A_init[2]),
                     R[4,3] / (1 - A_init[4] * A_init[3]), 
                     R[4,4] / (1 - A_init[4] * A_init[4])), 
                   ncol = 4, byrow = T)
    
    for(k in 1:n_i) {
        
        x_i_temp = matrix(c(x_i[k,]), ncol = 1)
        X_i[[k]] = diag(4) %x% x_i[k,]
        
        if(k == 1)  {
            D_i_temp = matrix(c(1, sum(before_b_i==2), sum(before_b_i==3), 
                                sum(before_b_i==4), sum(before_b_i==5)), 
                              nrow = 1, ncol = 5)
            D_i[[k]] = diag(4) %x% D_i_temp
            
            mean_vecY_i_k = D_i[[k]] %*% matrix(vec_alpha_i,ncol=1) + 
                X_i[[k]] %*% matrix(beta,ncol=1) + 
                D_i_omega[[k]] %*% matrix(omega_tilde,ncol=1)
            
            Y_i[k,] = rmvnorm(n=1, mean = mean_vecY_i_k, sigma = Gamma)
            
        } else {
            D_i_temp = matrix(c(1, sum(before_b_i==2) + sum(b_i[2:k]==2), 
                                sum(before_b_i==3) + sum(b_i[2:k]==3), 
                                sum(before_b_i==4) + sum(b_i[2:k]==4), 
                                sum(before_b_i==5) + sum(b_i[2:k]==5)), 
                              nrow = 1, ncol = 5)
            D_i[[k]] = diag(4) %x% D_i_temp
            
            A_curr = A_total[,b_i[k]]
            A_1 = diag(A_curr)
            
            nu_k = D_i[[k]] %*% matrix(vec_alpha_i,ncol=1) + 
                X_i[[k]] %*% matrix(beta,ncol=1) + 
                D_i_omega[[k]] %*% matrix(omega_tilde,ncol=1)
            nu_k_1 = D_i[[k-1]] %*% matrix(vec_alpha_i,ncol=1) + 
                X_i[[k-1]] %*% matrix(beta,ncol=1) + 
                D_i_omega[[k-1]] %*% matrix(omega_tilde,ncol=1)
            diff_vec = c(Y_i[k-1,] - nu_k_1)
            
            mean_vecY_i_k = nu_k + A_1 %*% matrix(diff_vec,ncol=1)
            
            Y_i[k,] = rmvnorm(n=1, mean = mean_vecY_i_k, sigma = R)
        }
    }
    # ----------------------------------------------------------------------
    
    # Verifying the rules --------------------------------------------------
    rules = Y[ Y[,'EID']==as.numeric(id_num), c('RBC_rule'), drop=F]
    rules = cbind(rules, 0)
    
    if(rbc_rule) {
        rbc_bleed_correct = c(rbc_bleed_correct, -1)
        
        if(2 %in% b_i) {
            first_bleed_ind = which(bleed_ind_i == 1)
            sim_bleed_ind = which(b_i == 2)
            # Check if bleed occurred before the RBC times (like we'd expect)
            if(2 %in% b_i[1:first_bleed_ind]){
                correct_bleed = T
                rbc_bleed_correct[length(rbc_bleed_correct)] = 1
            } 
        } 
    }
    
    if((1 %in% rules[,1]) & !(correct_bleed)) {
        rules[,1] = 0
    }
    
    data_format = rbind( data_format, cbind( id_num, Y_i, b_i, 
                                             z_i[,2],
                                             x_i[,1], 
                                             rules))
    initial_state_vec = c(initial_state_vec, b_i[1])
}

data_format = matrix(as.numeric(data_format), ncol = ncol(data_format))
colnames(data_format) = c('EID', 'hemo', 'hr', 'map','lactate', 
                          'b_true', 'RBC_ordered', 'n_RBC_admin', 
                          'RBC_rule', 'clinic_rule')

# Adding clinic rule -------------------------------------------------------
all_poss_combo = matrix(c(0,-1,
                          0, 1,
                          1, 1), nrow = 3, byrow = T)
subjects_2_3 = unique(data_format[data_format[,'b_true'] %in% c(2,3), "EID"])
subjects_2 = unique(data_format[data_format[,'b_true'] == 2, "EID"])

for(apc in 1:nrow(all_poss_combo)) {
    if(all_poss_combo[apc, 2] == -1) {
        # Clinic rule = -1
        num_rules = sum(clinic_rbc_combos[,3] == all_poss_combo[apc, 2])
        poss_EIDs = unique(data_format[data_format[,'RBC_rule'] == 0 & 
                                           !(data_format[,'EID'] %in% subjects_2_3), 'EID'])
        if(num_rules > length(poss_EIDs)) {
            clinic_assign = poss_EIDs
        } else {
            clinic_assign = sample(poss_EIDs, num_rules, replace = F)    
        }
        data_format[data_format[,"EID"] %in% clinic_assign, "clinic_rule"] = -1
    } else {
        # Clinic rule = 1
        num_rules = sum(clinic_rbc_combos[,2] == all_poss_combo[apc, 1] &
                            clinic_rbc_combos[,3] == all_poss_combo[apc, 2])
        poss_EIDs = unique(data_format[data_format[,"RBC_rule"] == all_poss_combo[apc, 1] & 
                                           data_format[,"EID"] %in% subjects_2, "EID"])
        clinic_assign = sample(poss_EIDs, num_rules, replace = F)
        
        data_format[data_format[,"EID"] %in% clinic_assign, "clinic_rule"] = 1
    }
}

# Impose missingness like real data ----------------------------------------
true_vitals = data_format[,c("hemo", "hr", "map", "lactate")]
colnames(true_vitals) = c('hm_true', 'hr_true', 'mp_true', 'la_true')
# data_format[!otype[,"hemo"], "hemo"] = NA
# data_format[!otype[,"hr"], "hr"] = NA
# data_format[!otype[,"map"], "map"] = NA
# data_format[!otype[,"lactate"], "lactate"] = NA
data_format = cbind(data_format, true_vitals)

# Print transition frequencies for the simulated data set ------------------
nTrans_sim = matrix(0, nrow = n_state, ncol = n_state)
for(i in unique(data_format[,"EID"])){
    subject <- data_format[data_format[,"EID"]==i,,drop=FALSE]
    for(k in 2:nrow(subject)) {
        nTrans_sim[subject[k-1, "b_true"], subject[k, "b_true"]] = 
            nTrans_sim[subject[k-1, "b_true"], subject[k, "b_true"]] + 1 
    }
}

change_states = c(nTrans_sim[1,2], nTrans_sim[1,4], nTrans_sim[2,3], nTrans_sim[2,4], 
                  nTrans_sim[3,1], nTrans_sim[3,2], nTrans_sim[3,4], nTrans_sim[4,2],
                  nTrans_sim[4,5], nTrans_sim[5,1], nTrans_sim[5,2], nTrans_sim[5,4])
print("All observed transitions: ")
print(nTrans_sim)

# Print summaries ----------------------------------------------------------
print(paste0(sum(rbc_bleed_correct == 1), " of ", length(rbc_bleed_correct), 
             " RBC rules were correct with the bleed event"))

print(paste0("Clinic Rule = (1, -1): ", "(", 
             length(unique(data_format[data_format[,"clinic_rule"] == 1, "EID"])),
             ", ", length(unique(data_format[data_format[,"clinic_rule"] == -1, "EID"])),
             ")"))

print("Initial state distribution")
print(table(initial_state_vec))
cat('\n')

# Save ---------------------------------------------------------------------
save(data_format, file=paste0('Data/sim_data_', df_num, '.rda'))

save(alpha_i_mat, file = paste0('Data/alpha_i_mat_', df_num, '.rda'))

save(bleed_indicator, file = paste0('Data/bleed_indicator_sim_', df_num, '.rda'))

if(df_num == 1) {
    save(Dn_omega_sim, file = paste0('Data/Dn_omega_sim.rda'))
}
# }

# # Visualize the noise --------------------------------------------------------
# EIDs = unique(data_format[,'EID'])
# simulation = T
# 
# makeTransparent = function(..., alpha=0.35) {
# 
#     if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
# 
#     alpha = floor(255*alpha)
#     newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
# 
#     .makeTransparent = function(col, alpha) {
#         rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
#     }
# 
#     newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
# 
#     return(newColor)
# 
# }
# 
# # New patients -----------------------------------------------------------------
# load('../Data_cleaning/Data/hr_map_names.rda')
# 
# pdf('Plots/initial_charts.pdf')
# panels = c(3, 1)
# inset_dim = c(0,-.18)
# par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
# for(i in EIDs[1:100]){
#     # print(which(EIDs == i))
#     indices_i = (data_format[,'EID']==i)
#     n_i = sum(indices_i)
#     t_grid = seq( 0, n_i, by=5)[-1]
#     rbc_times = which(data_format[indices_i, 'RBC_ordered'] != 0)
#     rbc_admin_times = which(diff(data_format[indices_i, 'n_RBC_admin']) > 0) + 1
# 
#     if(simulation) {
#         # Put this on the correct scale as the t_grid
#         b_i = data_format[ indices_i,'b_true']
#         to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
#         to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
#         to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
#         to_s4 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==4]
#         to_s5 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==5]
# 
#         if(b_i[1] == 1) {
#             to_s1 = c(to_s1, 1)
#         } else if(b_i[1] == 2) {
#             to_s2 = c(to_s2, 1)
#         } else if(b_i[1] == 3){
#             to_s3 = c(to_s3, 1)
#         } else if(b_i[1] == 4) {
#             to_s4 = c(to_s4, 1)
#         } else {
#             to_s5 = c(to_s5, 1)
#         }
# 
#         if(length(unique(b_i)) > 1) {
#             if(length(to_s1) > 0) {
#                 rect_coords = data.frame(s = 1, t = to_s1)
#             }
# 
#             if(length(to_s2) > 0) {
#                 s2_coords = data.frame(s = 2, t = to_s2)
#                 if(length(to_s1) > 0) {
#                     rect_coords = rbind(rect_coords, s2_coords)
#                 } else {
#                     rect_coords = s2_coords
#                 }
#             }
# 
#             if(length(to_s3) > 0) {
#                 s3_coords = data.frame(s = 3, t = to_s3)
#                 if(length(to_s1) > 0 || length(to_s2) > 0) {
#                     rect_coords = rbind(rect_coords, s3_coords)
#                 } else {
#                     rect_coords = s3_coords
#                 }
#             }
# 
#             if(length(to_s4) > 0) {
#                 s4_coords = data.frame(s = 4, t = to_s4)
#                 if(length(to_s1) > 0 || length(to_s2) > 0 || length(to_s3) > 0) {
#                     rect_coords = rbind(rect_coords, s4_coords)
#                 } else {
#                     rect_coords = s4_coords
#                 }
#             }
# 
#             if(length(to_s5) > 0) {
#                 s5_coords = data.frame(s = 5, t = to_s5)
#                 if(length(to_s1) > 0 || length(to_s2) > 0 || length(to_s3) > 0 || length(to_s4) > 0) {
#                     rect_coords = rbind(rect_coords, s5_coords)
#                 } else {
#                     rect_coords = s5_coords
#                 }
#             }
# 
#             if(!(n_i %in% rect_coords$t)) rect_coords = rbind(rect_coords, c(b_i[n_i], n_i))
#             # Add one row for visuals
#             rect_coords = rbind(rect_coords, c(b_i[n_i], n_i+1))
#             rect_coords$t = rect_coords$t - 1
#             rect_coords = rect_coords[order(rect_coords$t), ]
#             col_vec = c('dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray')[rect_coords$s]
#             col_vec = makeTransparent(col_vec, alpha = 0.35)
#         } else {
#             rect_coords = data.frame(s = rep(b_i[1], 2), t = c(1,n_i+1))
#             rect_coords$t = rect_coords$t - 1
#             rect_coords = rect_coords[order(rect_coords$t), ]
#             col_vec = c('dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray')[rect_coords$s]
#             col_vec = makeTransparent(col_vec, alpha = 0.35)
#         }
#     }
#     # HEART RATE & MAP ---------------------------------------------------------
#     hr_map_ylim = c(min(c(data_format[indices_i, 'hr'], data_format[indices_i, 'map'])),
#                     max(c(data_format[indices_i, 'hr'], data_format[indices_i, 'map'])))
# 
#     plot(1:n_i, data_format[indices_i, 'hr'], xlab='time', ylab=NA, ylim = hr_map_ylim,
#          main=paste0('HR and MAP: ', i, ', RBC Rule = ', mean(data_format[indices_i, 'RBC_rule'])),
#          col.main='green', col.axis='green', pch=20, cex=1, xaxt = 'n')
#     points(1:n_i, data_format[indices_i, 'map'], col = 'orange')
#     rect(xleft = rect_coords$t[-nrow(rect_coords)]-0.5,
#          ybottom = hr_map_ylim[1],
#          xright = rect_coords$t[-1]-0.5,
#          ytop = hr_map_ylim[2],
#          col = col_vec[-nrow(rect_coords)],
#          border = NA)
#     grid( nx=20, NULL, col='white')
#     abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
#     abline(v = rbc_admin_times, col = 'grey', lwd = 1)
#     xtick<-0:n_i
#     axis(side=1, at=xtick, labels = FALSE, col = 'green')
#     axis(side=1, at=t_grid, labels = t_grid, col = 'pink', col.axis = 'green')
# 
#     # HEMO & Lactate -----------------------------------------------------------
#     hemo_lact_ylim = c(min(c(data_format[indices_i, 'hemo'], data_format[indices_i, 'lactate'])),
#                        max(c(data_format[indices_i, 'hemo'], data_format[indices_i, 'lactate'])))
# 
#     plot(1:n_i, data_format[indices_i, 'hemo'], xlab='time', ylab=NA, ylim = hemo_lact_ylim,
#          main=paste0('hemo: ', i), col.main='green', col.axis='green', pch=20, cex=1)
#     points(1:n_i, data_format[indices_i, 'lactate'], col = 'orange')
#     rect(xleft = rect_coords$t[-nrow(rect_coords)]-0.5,
#          ybottom = hemo_lact_ylim[1],
#          xright = rect_coords$t[-1]-0.5,
#          ytop = hemo_lact_ylim[2],
#          col = col_vec[-nrow(rect_coords)],
#          border = NA)
#     grid( nx=20, NULL, col='white')
#     abline(v = rbc_times, col = 'darkorchid2', lwd = 1)
#     abline(v = rbc_admin_times, col = 'grey', lwd = 1)
# 
#     # Medication ---------------------------------------------------------------
#     med_i = Dn_omega_sim[[which(EIDs == i)]]
#     omega_i = omega_i_mat[[which(EIDs == i)]]
#     med_i_mat = stacked_chains = do.call( rbind, med_i)
# 
#     hr_med_i_mat = med_i_mat[seq(2, nrow(med_i_mat), by = 4), ]
#     map_med_i_mat = med_i_mat[seq(3, nrow(med_i_mat), by = 4), ]
# 
#     hr_mean_effect = hr_med_i_mat %*% omega_i
#     map_mean_effect = map_med_i_mat %*% omega_i
# 
#     hr_med_i_mat = hr_med_i_mat[, hr_map_names %in% c('hr_cont', 'hr_disc')]
#     map_med_i_mat = map_med_i_mat[, hr_map_names %in% c('map_cont', 'map_disc')]
# 
#     upp_down = c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
#                  -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
#                  -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
#                  -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
# 
#     hr_upp_down = upp_down[hr_map_names %in% c('hr_cont', 'hr_disc')]
#     map_upp_down = upp_down[hr_map_names %in% c('map_cont', 'map_disc')]
# 
#     hr_upp_i   = hr_med_i_mat[,hr_upp_down == 1]
#     hr_down_i  = hr_med_i_mat[,hr_upp_down == -1]
#     map_upp_i  = map_med_i_mat[,map_upp_down == 1]
#     map_down_i = map_med_i_mat[,map_upp_down == -1]
# 
#     total_hr_up  = rowSums(hr_upp_i)
#     total_hr_dn  = rowSums(hr_down_i)
#     total_map_up = rowSums(map_upp_i)
#     total_map_dn = rowSums(map_down_i)
# 
#     hr_map_ylim = c(min(total_hr_up, total_hr_dn, total_map_up, total_map_dn,
#                         hr_mean_effect, map_mean_effect),
#                     max(total_hr_up, total_hr_dn, total_map_up, total_map_dn,
#                         hr_mean_effect, map_mean_effect))
#     if(hr_map_ylim[1] == hr_map_ylim[2]) hr_map_ylim = c(0,1)
# 
#     plot(NULL, xlim=c(0.5,n_i+0.5), ylim=hr_map_ylim, main='Med. admin',
#          xlab='time', ylab=NA, col.main='green',
#          col.axis='green')
# 
#     if(simulation) {
#         rect(xleft = rect_coords$t[-nrow(rect_coords)],
#              ybottom = hr_map_ylim[1],
#              xright = rect_coords$t[-1],
#              ytop = hr_map_ylim[2],
#              col = col_vec[-nrow(rect_coords)],
#              border = NA)
#     }
#     points(x = 1:n_i, y = hr_mean_effect, xlab='time', ylab=NA,
#            col.main='green', col.axis='green',
#            col = 'aquamarine', pch = 16)
#     points(x = 1:n_i, y = map_mean_effect, xlab='time', ylab=NA,
#            col = 'orange', pch = 16)
#     lines(x = 1:n_i, y = hr_mean_effect, xlab='time', ylab=NA,
#           lwd=2, lty = 1, col = 'aquamarine')
#     lines(x = 1:n_i, y = map_mean_effect, xlab='time', ylab=NA,
#           lwd=2, lty = 1, col = 'orange')
# 
#     lines(x = 1:n_i, y = total_hr_up, xlab='time', ylab=NA,
#           lwd=1, lty = 2, col = 'aquamarine4')
#     lines(x = 1:n_i, y = total_map_up, xlab='time', ylab=NA,
#           lwd=1, lty = 3, col = 'darkolivegreen2')
#     lines(x = 1:n_i, y = total_hr_dn, xlab='time', ylab=NA,
#           lwd=1, lty = 4, col = 'deeppink')
#     lines(x = 1:n_i, y = total_map_dn, xlab='time', ylab=NA,
#           lwd=1, lty = 5, col = 'palevioletred')
#     legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
#             legend=c( 'HR effect', 'MAP effect'), pch=15, pt.cex=1.5,
#             col=c( 'aquamarine', 'orange'))
#     
# 
# }
# dev.off()