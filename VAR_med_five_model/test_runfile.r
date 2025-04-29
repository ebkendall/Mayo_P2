# ------------------------------------------------------------------------------
# Current values ---------------------------------------------------------------
# ------------------------------------------------------------------------------
alpha_tilde = c(  9, -1,  1, 0, 0,
                 85,  5, -5, 0, 0,
                 75, -5,  5, 0, 0,
                  5,  1, -1, 0, 0)

omega_tilde = 2*c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
                  -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
                  -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
                  -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)

# upsilon_alpha = diag(c(  4, 0.01, 0.01, 0.25, 0.25,
#                        100,    1,    1,   25,   25,
#                        100,    1,    1,   25,   25,
#                          1, 0.01, 0.01, 0.25, 0.25))
upsilon_alpha = diag(c(  4, 0.25, 0.25,  25,  25,
                       100,    9,    9, 100, 100,
                       100,    9,    9, 100, 100,
                         4, 0.25, 0.25,  25,  25))

upsilon_omega = diag(exp(rep(1, 84)))

# ------------------------------------------------------------------------------
# Priors (alpha_tilde, omega_tilde, upsilon_alpha) -----------------------------
# ------------------------------------------------------------------------------
alpha_0 = c( 9, -1,  1, 0, 0, 85,  5, -5, 0, 0,
            75, -5,  5, 0, 0,  5,  1, -1, 0, 0)
# sigma_alpha = diag(rep(20, length(alpha_0)))
sigma_alpha = diag(c(  4, 0.25, 0.25,  25,  25,
                     100,    9,    9, 100, 100,
                     100,    9,    9, 100, 100,
                       4, 0.25, 0.25,  25,  25))

omega_0 = 2*c(-1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
              -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
              -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
              -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1)
sigma_omega = diag(rep(20, length(omega_0)))

# nu_upsilon = 22
# psi_upsilon = diag(c(  4, 0.01, 0.01, 0.25, 0.25,
#                      100,    1,    1,   25,   25,
#                      100,    1,    1,   25,   25,
#                        1, 0.01, 0.01, 0.25, 0.25))
nu_upsilon = 50
psi_upsilon = diag(c(  4, 0.25, 0.25,  25,  25,
                     100,    9,    9, 100, 100,
                     100,    9,    9, 100, 100,
                       4, 0.25, 0.25,  25,  25))
psi_upsilon = (nu_upsilon - 20 - 1) * psi_upsilon

# ------------------------------------------------------------------------------
# Plotting samples from the priors ---------------------------------------------
# ------------------------------------------------------------------------------
library(mvtnorm, quietly = T)
library(LaplacesDemon, quietly = T)
library(ggplot2, quietly = T)
library(gridExtra, quietly = T)
set.seed(2025)

N = 2000
vital = c('Hemoglobin', 'Heart Rate', 'MAP', 'Lactate')
state_label = c(rep('S2', N), rep('S3', N), rep('S4', N), rep('S5', N))
color_choices = c("#E69F00", "#009E73", "#56B4E9", "#CC79A7")

pdf('Plots/re_prior_visual.pdf')
# alpha_tilde  -----------------------------------------------------------------
prior_sample_alpha_tilde = rmvnorm(N, mean = alpha_0, sigma = sigma_alpha)
plot_alpha_tilde_base = list()
plot_alpha_tilde_slope = list()
for(v in 1:length(vital)) {
    # baseline
    p_b = data.frame(state = rep('base', N),
                     samps = c(prior_sample_alpha_tilde[,5*(v-1)+1]))
    plot_alpha_tilde_base[[v]] = ggplot(p_b, aes(samps, fill = state)) + 
                                    geom_density(alpha = 0.2) + 
                                    labs(title = paste0(vital[v], ' Baseline Means'), 
                                         x = "", y = "") +
                                    theme(legend.position="none",
                                          plot.title = element_text(size=10))
    
    # slopes
    col_nums = (5*(v-1)+2):(5*v)
    p_s = data.frame(state = state_label,
                     samps = c(prior_sample_alpha_tilde[,col_nums]))
    
    plot_alpha_tilde_slope[[v]] = ggplot(p_s, aes(samps, fill = state)) + 
                                    geom_density(alpha = 0.2) + 
                                    labs(title = paste0(vital[v], ' Slope Means'),
                                         x = "", y = "") +
                                    theme(plot.title = element_text(size=10)) +
                                    scale_fill_manual(values = color_choices)
}

grid.arrange(plot_alpha_tilde_base[[1]], plot_alpha_tilde_slope[[1]], 
             plot_alpha_tilde_base[[2]], plot_alpha_tilde_slope[[2]],
             plot_alpha_tilde_base[[3]], plot_alpha_tilde_slope[[3]],
             plot_alpha_tilde_base[[4]], plot_alpha_tilde_slope[[4]], 
             nrow = 4, ncol = 2)

# upsilon_alpha ----------------------------------------------------------------
prior_sample_upsilon_alpha = matrix(nrow = N, ncol = ncol(psi_upsilon))
for(i in 1:N) {
    prior_sample_upsilon_alpha[i, ] = diag(rinvwishart(nu = nu_upsilon, S = psi_upsilon))
}

plot_upsilon_alpha_base = list()
plot_upsilon_alpha_slope = list()
for(v in 1:length(vital)) {
    # baseline
    p_b = data.frame(state = rep('base', N),
                     samps = c(prior_sample_upsilon_alpha[,5*(v-1)+1]))
    mean_b = round(mean(p_b$samps), digits = 4)
    med_b = round(median(p_b$samps), digits = 4)
    plot_upsilon_alpha_base[[v]] = ggplot(p_b, aes(samps, fill = state)) + 
        geom_density(alpha = 0.2) + 
        labs(title = paste0(vital[v], ' Baseline Variance'), 
             x = paste0("mean = ", mean_b, ", median = ", med_b), y = "") +
        theme(legend.position="none",
              plot.title = element_text(size=10),
              axis.title=element_text(size=8))
    
    col_nums = (5*(v-1)+2):(5*v)
    p_s = data.frame(state = state_label,
                     samps = c(prior_sample_upsilon_alpha[,col_nums]))
    
    plot_upsilon_alpha_slope[[v]] = list()
    for(j in 1:4) {
        p_s_j = p_s[p_s$state == unique(state_label)[j], ]
        mean_j = round(mean(p_s_j$samps), digits = 4)
        med_j = round(median(p_s_j$samps), digits = 4)
        plot_upsilon_alpha_slope[[v]][[j]] = ggplot(p_s_j, aes(samps, fill = state)) + 
            geom_density(alpha = 0.2) + 
            labs(title = paste0(vital[v], ' Slope Variance ', unique(state_label)[j]),
                 x = paste0("mean = ", mean_j, ", median = ", med_j), y = "") +
            theme(legend.position="none",
                  plot.title = element_text(size=10),
                  axis.title=element_text(size=8)) +
            scale_fill_manual(values = color_choices[j])
    }
    
    grid.arrange(plot_upsilon_alpha_base[[1]], 
                 plot_upsilon_alpha_slope[[v]][[1]], 
                 plot_upsilon_alpha_slope[[v]][[2]], 
                 plot_upsilon_alpha_slope[[v]][[3]], 
                 plot_upsilon_alpha_slope[[v]][[4]], 
                 nrow = 3, ncol = 2)
}

# grid.arrange(plot_upsilon_alpha_base[[1]], plot_upsilon_alpha_slope[[1]], 
#              plot_upsilon_alpha_base[[2]], plot_upsilon_alpha_slope[[2]],
#              plot_upsilon_alpha_base[[3]], plot_upsilon_alpha_slope[[3]],
#              plot_upsilon_alpha_base[[4]], plot_upsilon_alpha_slope[[4]], 
#              nrow = 4, ncol = 2)
# diag_variances = round(2 * diag(psi_upsilon)^2 / ((nu_upsilon - 20 - 1)^2 * (nu_upsilon - 20 - 3)), digits = 2)

# sampled alpha_i --------------------------------------------------------------
prior_sample_alpha_i = rmvnorm(N, mean = alpha_tilde, sigma = upsilon_alpha)
plot_alpha_i_base = list()
plot_alpha_i_slope = list()
for(v in 1:length(vital)) {
    # baseline
    p_b = data.frame(state = rep('base', N),
                     samps = c(prior_sample_alpha_i[,5*(v-1)+1]))
    plot_alpha_i_base[[v]] = ggplot(p_b, aes(samps, fill = state)) + 
        geom_density(alpha = 0.2) + 
        labs(title = paste0(vital[v], ' Baseline Random Effects'), 
             x = "", y = "") +
        theme(legend.position="none",
              plot.title = element_text(size=10))
    
    # slopes
    col_nums = (5*(v-1)+2):(5*v)
    p_s = data.frame(state = state_label,
                     samps = c(prior_sample_alpha_i[,col_nums]))
    
    plot_alpha_i_slope[[v]] = ggplot(p_s, aes(samps, fill = state)) + 
        geom_density(alpha = 0.2) + 
        labs(title = paste0(vital[v], ' Slope Random Effects'),
             x = "", y = "") +
        theme(plot.title = element_text(size=10)) +
        scale_fill_manual(values = color_choices)
}

grid.arrange(plot_alpha_i_base[[1]], plot_alpha_i_slope[[1]], 
             plot_alpha_i_base[[2]], plot_alpha_i_slope[[2]],
             plot_alpha_i_base[[3]], plot_alpha_i_slope[[3]],
             plot_alpha_i_base[[4]], plot_alpha_i_slope[[4]], 
             nrow = 4, ncol = 2)
dev.off()

