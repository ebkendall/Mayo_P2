source('mcmc_routine.r')

args = commandArgs(TRUE)
seed_num = as.numeric(args[1])
sampling_num = NULL

if(seed_num <= 3) {
    sampling_num = 1
} else if(seed_num > 3 & seed_num <= 6) {
    seed_num = seed_num - 3
    sampling_num = 2
} else if(seed_num > 6 & seed_num <= 9) {
    seed_num = seed_num - 6
    sampling_num = 3
} else {
    seed_num = seed_num - 9
    sampling_num = 4
}

set.seed(seed_num)
ind = seed_num

simulation = T
data_format = NULL

if(simulation) {
    steps  = 20000
    burnin =  5000

    trialNum = 4
    max_ind = 5
    sim_dat_num = 5
    
    load(paste0('Data_sim/use_data1_', sim_dat_num, '.rda'))
    data_format = use_data

    print(paste0('SIM: seed ', seed_num, ' samp ', 
                 sampling_num, ' trial ', trialNum))
} else {
    steps  = 50000
    burnin = 5000
    
    trialNum = 1
    max_ind = 5
    if(max_ind > 5) burnin = 0

    load('Data_updates/data_format.rda')
    print(paste0('REAL: seed ', seed_num, ' samp ', 
                 sampling_num, ' trial ', trialNum))
}

Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate', 'RBC_rule', 'clinic_rule')] 
EIDs = as.character(unique(data_format[,'EID']))

# Covariate on the mean process
x = data_format[,c('n_RBC_admin'), drop=F]
p = ncol(x)

# Covariate on the state process
if(simulation) {
    z = cbind(1, data_format[,c('RBC_ordered'), drop=F])
    m = ncol(z)
} else {
    # Edit to exponentially decay

}

# Parameter initialization -----------------------------------------------------
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

par = rep(0, max(par_index$vec_upsilon_omega))

par[par_index$vec_beta] = c(0.5, -2, 2, -0.5)

par[par_index$vec_alpha_tilde] = c( 9.57729783, -1,  1, 0, 0,
                                   88.69780576,  5, -5, 0, 0,
                                   79.74903940, -5,  5, 0, 0,
                                     5.2113319,  1, -1, 0, 0)
par[par_index$vec_sigma_upsilon] = c(diag(c(  3, 0.01, 0.01, 0.25, 0.25, 
                                            100,    1,    1,    1,    1, 
                                            100,    1,    1,    1,    1, 
                                              3, 0.01, 0.01, 0.25, 0.25)))

par[par_index$vec_A] = c(rep(1.5, 4), rep(-1, 4), rep(0.1, 4), rep(0, 4), rep(0.1, 4))

par[par_index$vec_R] = c(diag(c(9, 81, 81, 9)))

#    transitions:              1->2,         1->4,         2->3,         2->4, 
#                              3->1,         3->2,         3->4,         4->2, 
#                              4->5,         5->1,         5->2,         5->4,
par[par_index$vec_zeta] = c(-4.7405, 4.5, -5.2152,   1, -3.6473,-0.5, -3.1475, -0.2, 
                            -6.4459,  -1, -3.9404,   2, -4.2151,   1, -4.1778,  2.5, 
                            -3.0523,   0, -6.4459,-0.2, -4.2404, 3.5, -4.2151,    1)

par[par_index$vec_init] = c(-1, 0, -0.5, 0.1)

par[par_index$omega_tilde]= c(-1, -1,  1, -1,  1,  1, -1, -1, -1,  1,  1,  1,  1,
                              -1,  1,  1,  1, -1,  1, -1,  1, -1, -1, -1, -1, -1,
                              -1, -1, -1, -1,  1, -1,  1, -1, -1,  1,  1, -1,  1,
                              -1, -1,  1,  1, -1, -1, -1, -1, -1, -1,  1, -1,  1,
                               1, -1,  1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,
                              -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1, -1,
                              -1, -1,  1, -1, -1, -1, -1, -1, -1,  1)

par[par_index$vec_upsilon_omega] = rep(-4, length(par_index$vec_upsilon_omega))

# -----------------------------------------------------------------------------

if(simulation) {
    load(paste0('Data_sim/true_pars_', sim_dat_num, '.rda'))
    load(paste0('Data_sim/alpha_i_mat_', sim_dat_num, '.rda'))
    load(paste0('Data_sim/omega_i_mat_', sim_dat_num, '.rda'))
    load(paste0('Data_sim/Dn_omega_sim_', sim_dat_num, '.rda'))
    load(paste0('Data_sim/bleed_indicator_sim_', sim_dat_num,'.rda'))

    par = true_pars
    Dn_omega = Dn_omega_sim
    
    b_chain = data_format[, "b_true"]
} else {
    # ----------------------------------------------------------------------
    if(max_ind > 5) {
        prev_file = paste0('Model_out/mcmc_out_interm_', ind, '_', trialNum, 
                           'it', max_ind-5, '_df', df_num, '.rda')
        load(prev_file)
        
        par_temp = mcmc_out_temp$chain[nrow(mcmc_out_temp$chain), ]
        rownames(par_temp) = NULL
        par = par_temp
        
        b_chain = c(mcmc_out_temp$B_chain[nrow(mcmc_out_temp$B_chain), ])
        
        rm(mcmc_out_temp)
    }
    # ----------------------------------------------------------------------
    
    load('Data_updates/Dn_omega1.rda')
    load('Data_updates/all_EIDs.rda')
    Dn_omega_big = Dn_omega
    eid_index = which(all_EIDs %in% EIDs)
    Dn_omega = vector(mode = 'list', length = length(eid_index))
    for(jjj in 1:length(eid_index)) {
        Dn_omega[[jjj]] = Dn_omega_big[[eid_index[jjj]]]
    }
    rm(Dn_omega_big)
    
    bleed_indicator = b_ind_fnc(data_format)
}
# -----------------------------------------------------------------------------
A = list()
W = list()
B = list()

for(i in EIDs){
    if(simulation) {
        A[[i]] = alpha_i_mat[[which(EIDs == i)]]
        W[[i]] = omega_i_mat[[which(EIDs == i)]]
    } else {
        A[[i]] = matrix(par[par_index$vec_alpha_tilde], ncol =1)
        W[[i]] = matrix(par[par_index$omega_tilde], ncol =1)
    }
    
    if(max_ind > 5) {
        b_temp = b_chain[data_format[,"EID"] == as.numeric(i)]   
    } else {
        init_logit = par[par_index$vec_init]
        init_logit = c(0, init_logit)
        P_i = exp(init_logit) / sum(exp(init_logit))
        zeta = matrix(par[par_index$vec_zeta], nrow = 2)
        colnames(zeta) = c('(1) 1->2', '(2) 1->4','(3) 2->3', '(4) 2->4', '(5) 3->1', 
                           '(6) 3->2', '(7) 3->4','(8) 4->2', '(9) 4->5', '(10) 5->1', 
                           '(11) 5->2', '(12) 5->4')
        z_temp = z[Y[,"EID"] == i, ]
        n_i = nrow(z_temp)
        
        b_temp = NULL
        for(k in 1:n_i){
            if(k==1){
                b_temp = sample(1:5, size=1, prob=P_i)
            } else{
                q1   = exp(z_temp[k,, drop=F] %*% zeta[,  1, drop=F]) 
                q2   = exp(z_temp[k,, drop=F] %*% zeta[,  2, drop=F])
                q3   = exp(z_temp[k,, drop=F] %*% zeta[,  3, drop=F])
                q4   = exp(z_temp[k,, drop=F] %*% zeta[,  4, drop=F])
                q5   = exp(z_temp[k,, drop=F] %*% zeta[,  5, drop=F]) 
                q6   = exp(z_temp[k,, drop=F] %*% zeta[,  6, drop=F])
                q7   = exp(z_temp[k,, drop=F] %*% zeta[,  7, drop=F])
                q8   = exp(z_temp[k,, drop=F] %*% zeta[,  8, drop=F])
                q9   = exp(z_temp[k,, drop=F] %*% zeta[,  9, drop=F]) 
                q10  = exp(z_temp[k,, drop=F] %*% zeta[,  10, drop=F])
                q11  = exp(z_temp[k,, drop=F] %*% zeta[,  11, drop=F])
                q12  = exp(z_temp[k,, drop=F] %*% zeta[,  12, drop=F])
                
                Q = matrix(c(  1,   q1,  0,  q2,  0,
                               0,    1, q3,  q4,  0,
                              q5,   q6,  1,  q7,  0,
                               0,   q8,  0,   1, q9,
                             q10,  q11,  0, q12,  1), ncol=5, byrow=T)
                
                P_i = Q / rowSums(Q)
                # Sample the latent state sequence
                b_temp = c( b_temp, sample(1:5, size=1, prob=P_i[tail(b_temp,1),]))
            }
            
        }
    }
    
    B[[i]] = matrix(b_temp, ncol = 1)
}
# -----------------------------------------------------------------------------

print("Starting values for the chain")
print("alpha_tilde")
print(round(par[par_index$vec_alpha_tilde], 3))

print("A")
vec_A_t_logit = par[par_index$vec_A]
vec_A_t = (exp(vec_A_t_logit) - 1) / (1 + exp(vec_A_t_logit)) 
mat_A_t = matrix(vec_A_t, nrow = 4)
print(mat_A_t)

print("R")
R_t = matrix(par[par_index$vec_R], ncol = 4)
print(R_t)

print("zeta")
zed = matrix(par[par_index$vec_zeta], nrow = 2)
colnames(zed) = c('(1) 1->2', '(2) 1->4','(3) 2->3', '(4) 2->4', '(5) 3->1', 
                   '(6) 3->2', '(7) 3->4','(8) 4->2', '(9) 4->5', '(10) 5->1', 
                   '(11) 5->2', '(12) 5->4')
print(zed)

vec_A1 = par[par_index$vec_A]
scale_A1 = (exp(vec_A1) - 1) / (1 + exp(vec_A1)) 

diag_gamma = c(R_t[1,1] / (1 - scale_A1[1]^2), R_t[2,2] / (1 - scale_A1[2]^2),
               R_t[3,3] / (1 - scale_A1[3]^2), R_t[4,4] / (1 - scale_A1[4]^2),
               R_t[1,1] / (1 - scale_A1[5]^2), R_t[2,2] / (1 - scale_A1[6]^2),
               R_t[3,3] / (1 - scale_A1[7]^2), R_t[4,4] / (1 - scale_A1[8]^2),
               R_t[1,1] / (1 - scale_A1[9]^2), R_t[2,2] / (1 - scale_A1[10]^2),
               R_t[3,3] / (1 - scale_A1[11]^2), R_t[4,4] / (1 - scale_A1[12]^2),
               R_t[1,1] / (1 - scale_A1[13]^2), R_t[2,2] / (1 - scale_A1[14]^2),
               R_t[3,3] / (1 - scale_A1[15]^2), R_t[4,4] / (1 - scale_A1[16]^2),
               R_t[1,1] / (1 - scale_A1[17]^2), R_t[2,2] / (1 - scale_A1[18]^2),
               R_t[3,3] / (1 - scale_A1[19]^2), R_t[4,4] / (1 - scale_A1[20]^2))
print(round(sqrt(diag_gamma), 3))

s_time = Sys.time()
mcmc_out = mcmc_routine( par, par_index, A, W, B, Y, x, z, steps, burnin, ind, 
                         trialNum, Dn_omega, simulation, bleed_indicator, 
                         max_ind, df_num, sampling_num)
e_time = Sys.time() - s_time; print(e_time)