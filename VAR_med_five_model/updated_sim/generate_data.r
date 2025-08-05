library(mvtnorm, quietly=T)

# Parameter initialization -----------------------------------------------------
par_index = list()
par_index$alpha_tilde = 1:20
par_index$upsilon = 21:420
par_index$A = 421:424
par_index$R = 425:440
par_index$zeta = 441:452
par_index$init = 453:456

par = rep(0, max(do.call('c', par_index)))
par[par_index$alpha_tilde] = c( 50, -5,  5, -2,  2,
                               100, 10,-10,  2, -2,
                               100,-10, 10,  2, -2,
                                50,  5, -5, -2,  2)
par[par_index$upsilon] = c(diag(c(25, 4, 4, 1, 1, 
                                  25, 4, 4, 1, 1, 
                                  25, 4, 4, 1, 1, 
                                  25, 4, 4, 1, 1)))
par[par_index$A] = rep(0, 4)
par[par_index$R] = c(diag(c(9, 9, 9, 9)))
#    transitions:          1->2,    1->4,    2->3,    2->4, 
#                          3->1,    3->2,    3->4,    4->2, 
#                          4->5,    5->1,    5->2,    5->4
par[par_index$zeta] = c(-3.7405, -4.2152, -2.6473, -2.1475, 
                        -3.4459, -2.9404, -3.2151, -3.1778, 
                        -2.0523, -3.4459, -3.2404, -3.2151)
par[par_index$init] = c(0,0,0,0)

n_state = 5

# Parameter initialization -----------------------------------------------------
alpha_tilde = matrix(par[par_index$alpha_tilde], nrow = n_state)
Upsilon = matrix(par[par_index$upsilon], ncol = 20)
vec_A = c(par[par_index$A])
A_mat_scale = exp(vec_A) / (1 + exp(vec_A)) # Support (0,1)

# columns: hemo, hr, map, lactate
R = matrix(par[par_index$R], ncol = 4)

zeta = matrix(par[par_index$zeta], nrow = 1)
colnames(zeta) = c('(1) 1->2', '(2) 1->4','(3) 2->3', '(4) 2->4', '(5) 3->1', 
                   '(6) 3->2', '(7) 3->4','(8) 4->2', '(9) 4->5', '(10) 5->1', 
                   '(11) 5->2', '(12) 5->4')

init_logit = par[par_index$init]
init_logit = c(0, init_logit)

true_pars = par
save(true_pars, file = 'Data/true_pars.rda')
# ------------------------------------------------------------------------------

load('Data/data_format_train.rda')
EIDs = unique(data_format[,"EID"])
N = length(EIDs)
Y = data_format[, c('EID','hemo', 'hr', 'map', 'lactate')] 
rm(data_format)

# for(df_num in 1:100) {
args = commandArgs(TRUE)
df_num = as.numeric(args[1])
    
    print(paste0("Data set number: ", df_num))
    set.seed(df_num)
    
    data_format = NULL
    initial_state_vec = NULL
    alpha_i_mat = vector(mode = "list", length = N)
    
    for(i in 1:N){
        
        id_num = EIDs[i]
        
        n_i = sum(Y[,'EID']==as.numeric(id_num))
        m_i = n_i + rpois(n = 1, lambda = 50)
        
        # Sample a long string of states ---------------------------------------
        big_b_i = rep(NA, m_i)
        big_b_i[1] = 1
        marginal_prob = matrix(nrow = m_i, ncol = n_state)
        marginal_prob[1, ] = c(1,0,0,0,0)
        for(k in 2:m_i) {
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
            marginal_prob[k, ] = marginal_prob[k-1, ,drop=F] %*% P_i
        }
        
        # # Select n_i of the states for the true latent state sequence
        # b_i = tail(big_b_i, n_i)
        # before_b_i = big_b_i[1:(m_i - n_i + 1)]
        b_i = head(big_b_i, n_i)
        before_b_i = 1
        
        # Generate data --------------------------------------------------------
        Y_i = matrix(nrow = n_i, ncol = 4)
        D_i = vector(mode = 'list', length = n_i)
        
        vec_alpha_i = rmvnorm( n=1, mean=c(alpha_tilde), sigma=Upsilon)
        alpha_i_mat[[i]] = matrix(vec_alpha_i, ncol = 1)
        
        Gamma = matrix(c(R[1,1] / (1 - A_mat_scale[1] * A_mat_scale[1]), 
                         R[1,2] / (1 - A_mat_scale[1] * A_mat_scale[2]),
                         R[1,3] / (1 - A_mat_scale[1] * A_mat_scale[3]), 
                         R[1,4] / (1 - A_mat_scale[1] * A_mat_scale[4]),
                         R[2,1] / (1 - A_mat_scale[2] * A_mat_scale[1]), 
                         R[2,2] / (1 - A_mat_scale[2] * A_mat_scale[2]),
                         R[2,3] / (1 - A_mat_scale[2] * A_mat_scale[3]), 
                         R[2,4] / (1 - A_mat_scale[2] * A_mat_scale[4]),
                         R[3,1] / (1 - A_mat_scale[3] * A_mat_scale[1]), 
                         R[3,2] / (1 - A_mat_scale[3] * A_mat_scale[2]),
                         R[3,3] / (1 - A_mat_scale[3] * A_mat_scale[3]), 
                         R[3,4] / (1 - A_mat_scale[3] * A_mat_scale[4]),
                         R[4,1] / (1 - A_mat_scale[4] * A_mat_scale[1]), 
                         R[4,2] / (1 - A_mat_scale[4] * A_mat_scale[2]),
                         R[4,3] / (1 - A_mat_scale[4] * A_mat_scale[3]), 
                         R[4,4] / (1 - A_mat_scale[4] * A_mat_scale[4])), 
                       ncol = 4, byrow = T)
        
        for(k in 1:n_i) {
            if(k == 1) {
                D_i_temp = matrix(c(1, sum(before_b_i==2), sum(before_b_i==3), 
                                    sum(before_b_i==4), sum(before_b_i==5)), 
                                  nrow = 1, ncol = 5)
                D_i[[k]] = diag(4) %x% D_i_temp
                
                mean_vecY_i_k = D_i[[k]] %*% matrix(vec_alpha_i,ncol=1) 
                
                Y_i[k,] = rmvnorm(n=1, mean = mean_vecY_i_k, sigma = Gamma)
                
            } else {
                D_i_temp = matrix(c(1, sum(before_b_i==2) + sum(b_i[2:k]==2), 
                                    sum(before_b_i==3) + sum(b_i[2:k]==3), 
                                    sum(before_b_i==4) + sum(b_i[2:k]==4), 
                                    sum(before_b_i==5) + sum(b_i[2:k]==5)), 
                                  nrow = 1, ncol = 5)
                D_i[[k]] = diag(4) %x% D_i_temp
                
                A_1 = diag(A_mat_scale)
                
                nu_k   = D_i[[k]] %*% matrix(vec_alpha_i,ncol=1) 
                nu_k_1 = D_i[[k-1]] %*% matrix(vec_alpha_i,ncol=1)
                diff_vec = c(Y_i[k-1,] - nu_k_1)
                
                mean_vecY_i_k = nu_k + A_1 %*% matrix(diff_vec,ncol=1)
                
                Y_i[k,] = rmvnorm(n=1, mean = mean_vecY_i_k, sigma = R)
            }
        }
        # ----------------------------------------------------------------------
        
        data_format = rbind( data_format, cbind( id_num, Y_i, b_i))
        initial_state_vec = c(initial_state_vec, b_i[1])
    }
    data_format = matrix(as.numeric(data_format), ncol = ncol(data_format))
    colnames(data_format) = c( 'EID', 'hemo', 'hr', 'map','lactate', 'b_true')
    
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
    
    print("Initial state distribution")
    print(table(initial_state_vec))
    cat('\n')
    
    # Save ---------------------------------------------------------------------
    save(data_format, file=paste0('Data/sim_data_', df_num, '.rda'))
    
    save(alpha_i_mat, file = paste0('Data/alpha_i_mat_', df_num, '.rda'))
# }

# # Visualize the noise --------------------------------------------------------
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
# pdf('Plots/initial_charts.pdf')
# panels = c(4, 1)
# inset_dim = c(0,-.18)
# par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
# for(i in EIDs[1:100]){
# 
#     indices_i = (data_format[,'EID']==i)
#     n_i = sum(indices_i)
#     t_grid = seq( 0, n_i, by=5)[-1]
# 
#     # Put this on the correct scale as the t_grid
#     b_i = data_format[ indices_i,'b_true']
#     to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
#     to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
#     to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
#     to_s4 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==4]
#     to_s5 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==5]
# 
#     if(b_i[1] == 1) {
#         to_s1 = c(to_s1, 1)
#     } else if(b_i[1] == 2) {
#         to_s2 = c(to_s2, 1)
#     } else if(b_i[1] == 3){
#         to_s3 = c(to_s3, 1)
#     } else if(b_i[1] == 4) {
#         to_s4 = c(to_s4, 1)
#     } else {
#         to_s5 = c(to_s5, 1)
#     }
# 
#     if(length(unique(b_i)) > 1) {
#         if(length(to_s1) > 0) {
#             rect_coords = data.frame(s = 1, t = to_s1)
#         }
# 
#         if(length(to_s2) > 0) {
#             s2_coords = data.frame(s = 2, t = to_s2)
#             if(length(to_s1) > 0) {
#                 rect_coords = rbind(rect_coords, s2_coords)
#             } else {
#                 rect_coords = s2_coords
#             }
#         }
# 
#         if(length(to_s3) > 0) {
#             s3_coords = data.frame(s = 3, t = to_s3)
#             if(length(to_s1) > 0 || length(to_s2) > 0) {
#                 rect_coords = rbind(rect_coords, s3_coords)
#             } else {
#                 rect_coords = s3_coords
#             }
#         }
# 
#         if(length(to_s4) > 0) {
#             s4_coords = data.frame(s = 4, t = to_s4)
#             if(length(to_s1) > 0 || length(to_s2) > 0 || length(to_s3) > 0) {
#                 rect_coords = rbind(rect_coords, s4_coords)
#             } else {
#                 rect_coords = s4_coords
#             }
#         }
# 
#         if(length(to_s5) > 0) {
#             s5_coords = data.frame(s = 5, t = to_s5)
#             if(length(to_s1) > 0 || length(to_s2) > 0 || length(to_s3) > 0 || length(to_s4) > 0) {
#                 rect_coords = rbind(rect_coords, s5_coords)
#             } else {
#                 rect_coords = s5_coords
#             }
#         }
# 
#         if(!(n_i %in% rect_coords$t)) rect_coords = rbind(rect_coords, c(b_i[n_i], n_i))
#         # Add one row for visuals
#         rect_coords = rbind(rect_coords, c(b_i[n_i], n_i+1))
#         rect_coords$t = rect_coords$t - 1
#         rect_coords = rect_coords[order(rect_coords$t), ]
#         col_vec = c('dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray')[rect_coords$s]
#         col_vec = makeTransparent(col_vec, alpha = 0.35)
#     } else {
#         rect_coords = data.frame(s = rep(b_i[1], 2), t = c(1,n_i+1))
#         rect_coords$t = rect_coords$t - 1
#         rect_coords = rect_coords[order(rect_coords$t), ]
#         col_vec = c('dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray')[rect_coords$s]
#         col_vec = makeTransparent(col_vec, alpha = 0.35)
#     }
# 
#     # HEART RATE & MAP ---------------------------------------------------------
#     hr_map_ylim = c(min(c(data_format[indices_i, 'hr'], data_format[indices_i, 'map'])),
#                     max(c(data_format[indices_i, 'hr'], data_format[indices_i, 'map'])))
# 
#     plot(1:n_i, data_format[indices_i, 'hr'], xlab='time', ylab=NA, ylim = hr_map_ylim,
#          main=paste0('HR and MAP: ', i),
#          col.main='green', col.axis='green', pch=20, cex=1, xaxt = 'n')
#     points(1:n_i, data_format[indices_i, 'map'], col = 'orange')
#     rect(xleft = rect_coords$t[-nrow(rect_coords)]-0.5,
#          ybottom = hr_map_ylim[1],
#          xright = rect_coords$t[-1]-0.5,
#          ytop = hr_map_ylim[2],
#          col = col_vec[-nrow(rect_coords)],
#          border = NA)
#     grid( nx=20, NULL, col='white')
#     xtick<-0:n_i
#     axis(side=1, at=xtick, labels = FALSE, col = 'green')
#     axis(side=1, at=t_grid, labels = t_grid, col = 'pink', col.axis = 'green')
# 
#     # HEMO & Lactate -----------------------------------------------------------
#     hemo_lact_ylim = c(min(c(data_format[indices_i, 'hemo'], data_format[indices_i, 'lactate'])),
#                        max(c(data_format[indices_i, 'hemo'], data_format[indices_i, 'lactate'])))
# 
#     plot(1:n_i, data_format[indices_i, 'hemo'], xlab='time', ylab=NA, ylim = hemo_lact_ylim,
#          main=paste0('HEMO and LACT: ', i), col.main='green', col.axis='green', pch=20, cex=1)
#     points(1:n_i, data_format[indices_i, 'lactate'], col = 'orange')
#     rect(xleft = rect_coords$t[-nrow(rect_coords)]-0.5,
#          ybottom = hemo_lact_ylim[1],
#          xright = rect_coords$t[-1]-0.5,
#          ytop = hemo_lact_ylim[2],
#          col = col_vec[-nrow(rect_coords)],
#          border = NA)
#     grid( nx=20, NULL, col='white')
# }
# dev.off()
