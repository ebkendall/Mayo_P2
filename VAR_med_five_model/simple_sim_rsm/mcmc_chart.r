library(matrixStats)
library(plotrix)

args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])
it_num = 3
states_per_step = 1
steps_per_it = 1
S = 2


seed_list = c(1:3)

# Mode of the state sequences -------------------------------------------------
Mode <- function(x) {
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
}

# Load the model output -------------------------------------------------------
B_chain = NULL
it_seq = 1:it_num

state_results = vector(mode = 'list', length = length(seed_list))

for(s in 1:length(seed_list)) {
    seed_num = seed_list[s]
    for(it in it_seq) {
        file_name = paste0('Model_out/mcmc_out_',seed_num,'_', 'it',
                           it, '_samp', sampling_num,
                           '_', states_per_step, '_', steps_per_it,'.rda')
        load(file_name)
        print(file_name)
        
        if(it == 1) {
            B_chain   = mcmc_out$B_chain[1000:2000, ]
        } else {
            B_chain   = rbind(B_chain, mcmc_out$B_chain)
        }
        rm(mcmc_out)
    }    
    
    state_results[[s]] = matrix(nrow = S+1, ncol = ncol(B_chain))
    for(jj in 1:S) {
        state_results[[s]][jj, ] = apply(B_chain, 2, function(x,jj){sum(x == jj)}, jj)
    }
    state_results[[s]][S+1, ] = apply(B_chain, 2, Mode)
    rm(B_chain)
}

combo_counts = state_results[[1]]
for(i in 2:length(state_results)) {
    combo_counts = combo_counts + state_results[[i]]
}

all_seeds_state_mode = matrix(nrow = length(seed_list)+1, ncol = ncol(combo_counts))
for(i in 1:length(seed_list)) {
    all_seeds_state_mode[i, ] = state_results[[i]][S+1,]
}
all_seeds_state_mode[length(seed_list)+1, ] = apply(combo_counts, 2, which.max)

load('Data/data_format1.rda')
y = data_format[,"y"]
ids = data_format[,"id"]
EIDs = unique(data_format[,"id"])
ss_truth = data_format[,"state"]
#  ----------------------------------------------------------------------------
# state_seq_mode = apply(B_chain, 2, Mode)
# ------------------------------------------------------------------------------
# Function to change transparency of colors # ----------------------------------
# ------------------------------------------------------------------------------
makeTransparent = function(..., alpha=0.35) {
    
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    
    alpha = floor(255*alpha)  
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    
    .makeTransparent = function(col, alpha) {
        rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    
    return(newColor)
    
}


print("Summary of identifying correct states with mode")
for(s in 1:(length(seed_list)+1)) {
    if(s <= length(seed_list)) {
        print(paste0("Seed ", s))
    } else {
        print("Combo")
    }
    
    if(length(ss_truth) != length(all_seeds_state_mode[s, ])) {
        print("ERROR")
    } else {
        print(sum(ss_truth == all_seeds_state_mode[s, ]) / length(all_seeds_state_mode[s, ]))   
    }
}

print("Average (across subject) of proportion of correct states with mode")
mode_correctness = matrix(nrow = length(EIDs), ncol = length(seed_list)+1)
for(s in 1:(length(seed_list)+1)) {
    prop_sub = rep(0, length(EIDs))
    for(j in 1:length(EIDs)) {
        sub_ind_j = which(data_format[,"id"] == EIDs[j])
        ss_truth_j = ss_truth[sub_ind_j]
        state_seq_mode_j = all_seeds_state_mode[s, sub_ind_j]
        
        prop_sub[j] = sum(ss_truth_j == state_seq_mode_j) / length(state_seq_mode_j)
    }

    if(s <= length(seed_list)) {
        print(paste0("Seed ", s))
    } else {
        print("Combo")
    }
    
    print(summary(prop_sub))
    mode_correctness[,s] = prop_sub
    
    eid_poor = NULL
    eid_poor = EIDs[prop_sub < 0.9]
    print("EIDs with < 90% of correct state identification")
    print(eid_poor)    
}

pdf_title = paste0('trace_plot_par_fix_samp', sampling_num, '_', states_per_step, 
                   '_', steps_per_it,'.pdf')
pdf(pdf_title)
par(mfrow=c(2,2))
for(i in 1:(length(seed_list)+1)) {
    plot_title = NULL
    if(i <= length(seed_list)) {
        plot_title = paste0("Seed ", i)
    } else {
        plot_title = "Combo" 
    }
    
    boxplot(mode_correctness[,i], main = plot_title, ylab = "percent accurate")
    abline(h = sum(ss_truth == all_seeds_state_mode[i, ]) / length(all_seeds_state_mode[i, ]), col = 'green')
}
dev.off()

# # ------------------------------------------------------------------------------ 
# # Model evaluation plots -------------------------------------------------------
# # ------------------------------------------------------------------------------
# pdf_title = paste0('Plots/chart_', seed_num, "_samp",  sampling_num, '_sim.pdf') 
# pdf(pdf_title)
# panel_dim = c(4,1)
# inset_dim = c(0,-.18)
# par(mfrow=panel_dim, mar=c(2,4,2,4), bg='black', fg='green')
# for(i in EIDs){
#     print(which(EIDs == i))
#     indices_i = (data_format[,'EID']==i)
#     n_i = sum(indices_i)
#     
#     t_grid = round(data_format[indices_i, 'time'] / 60, digits = 3)
#     t_grid_bar = 1:length(t_grid)
#     rbc_times_bar = which(data_format[data_format[,'EID']==i, 'RBC_ordered'] != 0)
#     if(simulation) {
#         rbc_admin = c(head(data_format[data_format[,'EID']==i, "n_RBC_admin"], 1),
#                       diff(data_format[data_format[,'EID']==i, "n_RBC_admin"]))
#         rbc_admin_times_bar = which(rbc_admin != 0)
#     } else {
#         rbc_admin_times_bar = which(data_format[data_format[,'EID']==i, 'RBC_admin'] != 0)   
#     }
#     rbc_times = t_grid[rbc_times_bar]
#     rbc_admin_times = t_grid[rbc_admin_times_bar]
#     
#     # Put this on the correct scale as the t_grid
#     b_i = data_format[ indices_i,'b_true']
#     to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
#     to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
#     
#     if(b_i[1] == 1) {
#         to_s1 = c(to_s1, 1)
#     } else if(b_i[1] == 2) {
#         to_s2 = c(to_s2, 1)
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
#             if(length(to_s1) > 0 || length(to_s2) > 0 || 
#                length(to_s3) > 0 || length(to_s4) > 0) {
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
#         col_vec = c('dodgerblue', 'firebrick1')[rect_coords$s]
#         col_vec = makeTransparent(col_vec, alpha = 0.35)   
#     } else {
#         rect_coords = data.frame(s = rep(b_i[1], 2), t = c(1,n_i+1))
#         rect_coords$t = rect_coords$t - 1
#         rect_coords = rect_coords[order(rect_coords$t), ]
#         col_vec = c('dodgerblue', 'firebrick1')[rect_coords$s]
#         col_vec = makeTransparent(col_vec, alpha = 0.35)  
#     }
#     
#     pb = barplot(rbind(colMeans(B_chain[, indices_i] == 1),
#                        colMeans(B_chain[, indices_i] == 2)), 
#                  col=c( 'dodgerblue', 'firebrick1'), 
#                  xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
#     
#     # Make a new plot to add the background color
#     plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
#          xlab='time', ylab=NA, xaxt='n', col.main='green',
#          col.axis='green')
#     
#     rect(xleft = rect_coords$t[-nrow(rect_coords)], 
#          ybottom = hr_map_ylim[1], 
#          xright = rect_coords$t[-1], 
#          ytop = hr_map_ylim[2],
#          col = col_vec[-nrow(rect_coords)],
#          border = NA)
#     
#     plotCI( x = pb, y=colMeans(Hr_chain[, indices_i, drop=F]), ui=hr_upper, li=hr_lower,
#             main=title_name,
#             xlab='time', ylab=NA, xaxt='n', col.main='green',
#             col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
#             xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim, add =T) 
#     plotCI( x = pb, y=colMeans(Map_chain[, indices_i, drop=F]), ui=bp_upper, li=bp_lower,
#             main=title_name,
#             xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
#             col = 'orange',
#             xlim = range(pb) + c(-0.5,0.5), add = T) 
#     legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
#             legend=c( 'HR', 'MAP'), pch=15, pt.cex=1.5, 
#             col=c( 'aquamarine', 'orange'))
#     grid( nx=20, NULL, col='white')
#     axis( side=1, at=pb, col.axis='green', labels=t_grid)
#     
#     abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
#     abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
#     
#     # BAR PLOTS --------------------------------------------------------------
#     barplot(rbind(colMeans(B_chain[, indices_i] == 1),
#                   colMeans(B_chain[, indices_i] == 2),
#                   colMeans(B_chain[, indices_i] == 3),
#                   colMeans(B_chain[, indices_i] == 4),
#                   colMeans(B_chain[, indices_i] == 5)), 
#             col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
#             xlab='time', space=0, col.main='green', border=NA,
#             xlim=range(pb) + c(-0.5,0.5)) 
#     grid( nx=20, NULL, col='white')
#     legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
#             legend=c( 'Baseline', 'State 2', 'State 3', 'State 4', 'State 5'), 
#             pch=15, pt.cex=1.5, 
#             col=c( 'dodgerblue', 'firebrick1', 'yellow2','green', 'darkgray'))
#     legend( 'topleft', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
#             legend=c( 'RBC order', 'RBC admin'), pch=15, pt.cex=1.5, 
#             col=c( 'darkorchid1', 'aquamarine'))				
#     axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
#     axis( side=2, at=0:1, col.axis='green')
#     
#     
#     abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
#     abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
# }
# dev.off()
