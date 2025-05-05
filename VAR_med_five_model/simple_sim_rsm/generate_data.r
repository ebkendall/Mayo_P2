library(mvtnorm, quietly=T)

# Parameter initialization -----------------------------------------------------
par_index = list()
par_index$alpha = 1:12
par_index$zeta = 13:16
par_index$diag_R = 17:20
par_index$init = 21:22

par = rep(0, tail(par_index$init, 1))
par[par_index$alpha] = c( 50, -3,  3,
                         100,  5, -5,
                         100, -5,  5,
                          50,  3, -3)
par[par_index$zeta] = c(-2, -1, -1.5, -1.5)
par[par_index$diag_R] = c(0, 0, 0, 0)
par[par_index$init] = c(0, 0)

N = 500
n_state = 3

# Defining parameter objects ---------------------------------------------------
alpha = matrix(par[par_index$alpha], nrow = 3)
zeta = par[par_index$zeta]
R = diag(exp(par[par_index$diag_R]))
logit_prob = c(1, exp(par[par_index$init]))
init_prob = logit_prob / sum(logit_prob)

qz = exp(zeta)
Q = matrix(c(    1, qz[1],     0,
                 0,     1, qz[2],
             qz[3], qz[4],     1), ncol=3, byrow=T)

P = Q / rowSums(Q)

# Simulate multiple datasets ---------------------------------------------------
for(seed_num in 1:100) {
    
    print(paste0("Data set number: ", seed_num))
    set.seed(seed_num)
    
    if(seed_num == 1) {
        pdf("Plots/sim_data.pdf")
        par(mfrow=c(2,1))
    }
    data_format = NULL
    for(i in 1:N) {
        n_i = rpois(n = 1, lambda = 100)
        b_i = rep(NA, n_i)
        y_i = matrix(nrow = n_i, ncol = 4)
        
        # Sample the latent states
        for(k in 1:n_i) {
            if(k == 1) {
                b_i[k] = sample(1:n_state, size = 1, prob = init_prob)
            } else {
                b_i[k] = sample(1:n_state, size = 1, prob=P[b_i[k-1],])
            }
        }

        # Sample the observations
        g_0 = rep(0, ncol(y_i))
        max_t = 20
        for(k in 1:n_i) {
            mean_k = rep(0, ncol(y_i))
            if(k == 1) {
                t_2 = 0
                t_3 = 0
                if(b_i[k] == 2) {
                    t_3 = sample(0:max_t, size = 1)
                    if(t_3 > 0) {
                        t_2 = sample(2:max_t, size = 1)
                    } else {
                        t_2 = sample(1:max_t, size = 1)
                    }
                } else if(b_i[k] == 3) {
                    t_2 = sample(1:max_t, size = 1)
                    t_3 = sample(1:max_t, size = 1)
                }
                g_0 = alpha[1,] + t_2 * alpha[2, ] + t_3 * alpha[3, ]
                mean_k = g_0
            } else {
                mean_k = g_0 + sum(b_i[2:k] == 2) * alpha[2, ] + sum(b_i[2:k] == 3) * alpha[3, ]
            }
            
            y_i[k, ] = rmvnorm(n = 1, mean = mean_k, sigma = R)
        }
        
        sub_i = cbind(i, y_i, b_i, t_2, t_3)
        data_format = rbind(data_format, sub_i)
        
        if(seed_num == 1 & i <= 50) {
            y_bounds = c(min(c(y_i)), max(c(y_i)))
            plot(1:n_i, y_i[,1], ylim = y_bounds, ylab = "y", 
                 xlab = paste0("Initial state = ", b_i[1], ", t2 = ", t_2, ", t3 = ", t_3))
            points(1:n_i, y_i[,2], col = 'red')
            points(1:n_i, y_i[,3], col = 'green')
            points(1:n_i, y_i[,4], col = 'purple')    
        }
    }
    if(seed_num == 1) dev.off()
    
    colnames(data_format) = c("id", "y1", "y2", "y3", "y4", "state", "t2", "t3")
    save(data_format, file = paste0('Data/data_format', seed_num, '.rda'))    

    # Transition frequencies for the simulated data set ------------------------
    nTrans_sim = matrix(0, nrow = n_state, ncol = n_state)
    for(i in unique(data_format[,"id"])){
    	subject <- data_format[data_format[,"id"]==i,,drop=FALSE]
    	for(k in 2:nrow(subject)) {
    	    nTrans_sim[subject[k-1, "state"], subject[k, "state"]] = 
    	        nTrans_sim[subject[k-1, "state"], subject[k, "state"]] + 1 
    	}
    }
    
    change_states = c(nTrans_sim[1,2], nTrans_sim[2,3], nTrans_sim[3,1], nTrans_sim[3,2])
    print("All observed transitions: ")
    print(nTrans_sim)
    cat('Transition fequencies = ', change_states / sum(change_states),'\n')
    cat('Transition counts     = ', change_states,'\n')
}
