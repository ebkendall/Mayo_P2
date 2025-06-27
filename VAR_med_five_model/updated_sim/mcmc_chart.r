index_seeds = c(1,3,4,5)
trialNum = 1
it_num = 1
S = 5

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
        file_name = paste0('Model_out/mcmc_out_',seed, '_it_', it,'.rda')
        file_name = paste0('Model_out/mcmc_out_',trialNum, '_', seed, 'it', it,'_sim.rda')
        if(file.exists(file_name)) {
            load(file_name)
            print(paste0(seed, ": ", file_name))
            
            par_index = mcmc_out$par_index
            
            if(it == 1) {
                ind = ind + 1
                B_chain = mcmc_out$B_chain[500:1000, ]
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
for(seed in index_seeds) {
    
    load(paste0('Data/sim_data_', seed, '.rda'))
    ss_true = data_format[,"b_true"]
    
    if(length(ss_true) != length(state_results[[seed]][S+1, ])) {
        print(paste0("ERROR (seed ", seed, ")"))
    } else {
        prop_all = sum(ss_true == state_results[[seed]][S+1, ]) / length(state_results[[seed]][S+1, ])
        print(paste0("Across all subjects (seed ", seed, ") = ", prop_all))
    }
    
    EIDs = unique(data_format[,"EID"])
    
    prop_sub = rep(0, length(EIDs))
    
    for(j in 1:length(EIDs)) {
        sub_ind_j = which(data_format[,"EID"] == EIDs[j])
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

# Plot Chart Plots for Seed 1 data ---------------------------------------------
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

seed_focus = 1
load(paste0('Data/sim_data_', seed_focus, '.rda'))
EIDs = unique(data_format[,"EID"])

pdf_title = paste0('Plots/chart_plot.pdf')
pdf(pdf_title)
panel_dim = c(3,1)
inset_dim = c(0,-.18)
par(mfrow=panel_dim, mar=c(2,4,2,4), bg='black', fg='green')
set.seed(2025)
eid_plot = sample(EIDs, size = 25)
for(i in eid_plot){
    
    indices_i = (data_format[,'EID']==i)
    n_i = sum(indices_i)
    
    t_grid = t_grid_bar = 1:n_i
    
    # Put this on the correct scale as the t_grid
    b_i = data_format[ indices_i,"b_true"]
    to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
    to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
    to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
    to_s4 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==4]
    to_s5 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==5]
    
    if(b_i[1] == 1) {
        to_s1 = c(to_s1, 1)
    } else if(b_i[1] == 2) {
        to_s2 = c(to_s2, 1)
    } else if(b_i[1] == 3) {
        to_s3 = c(to_s3, 1)
    } else if(b_i[1] == 4) {
        to_s4 = c(to_s4, 1)
    } else {
        to_s5 = c(to_s5, 1)
    }
    
    if(length(unique(b_i)) > 1) {
        if(length(to_s1) > 0) {
            rect_coords = data.frame(s = 1, t = to_s1)
        }
        
        if(length(to_s2) > 0) {
            s2_coords = data.frame(s = 2, t = to_s2)
            if(length(to_s1) > 0) {
                rect_coords = rbind(rect_coords, s2_coords)
            } else {
                rect_coords = s2_coords
            }
        }
        
        if(length(to_s3) > 0) {
            s3_coords = data.frame(s = 3, t = to_s3)
            if(length(to_s1) > 0 || length(to_s2) > 0) {
                rect_coords = rbind(rect_coords, s3_coords)
            } else {
                rect_coords = s3_coords
            }
        }
        
        if(length(to_s4) > 0) {
            s4_coords = data.frame(s = 4, t = to_s4)
            if(length(to_s1) > 0 || length(to_s2) > 0 || length(to_s3) > 0) {
                rect_coords = rbind(rect_coords, s4_coords)
            } else {
                rect_coords = s4_coords
            }
        }

        if(length(to_s5) > 0) {
            s5_coords = data.frame(s = 5, t = to_s5)
            if(length(to_s1) > 0 || length(to_s2) > 0 ||
               length(to_s3) > 0 || length(to_s4) > 0) {
                rect_coords = rbind(rect_coords, s5_coords)
            } else {
                rect_coords = s5_coords
            }
        }
        
        if(!(n_i %in% rect_coords$t)) rect_coords = rbind(rect_coords, c(b_i[n_i], n_i))
        # Add one row for visuals
        rect_coords = rbind(rect_coords, c(b_i[n_i], n_i+1))
        rect_coords$t = rect_coords$t - 1
        rect_coords = rect_coords[order(rect_coords$t), ]
        col_vec = c('dodgerblue', 'firebrick1', 'yellow2', 
                    'green', 'darkgray')[rect_coords$s]
        col_vec = makeTransparent(col_vec, alpha = 0.35)
    } else {
        rect_coords = data.frame(s = rep(b_i[1], 2), t = c(1,n_i+1))
        rect_coords$t = rect_coords$t - 1
        rect_coords = rect_coords[order(rect_coords$t), ]
        col_vec = c('dodgerblue', 'firebrick1', 'yellow2', 
                    'green', 'darkgray')[rect_coords$s]
        col_vec = makeTransparent(col_vec, alpha = 0.35)  
    }
    
    
    pb = barplot(state_results[[seed_focus]][1:S, indices_i], 
                 col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                 xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F)
    
    # TRUE STATE SEQ -----------------------------------------------------------
    plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=c(min(data_format[indices_i,c("hr", "map")]),
                                                    max(data_format[indices_i,c("hr", "map")])), 
         main=paste0("EID = ", i, " truth"), xlab='time', ylab=NA, xaxt='n', 
         col.main='green', col.axis='green')
    
    # Vitals 2 & 3 -------------------------------------------------------------
    rect(xleft = rect_coords$t[-nrow(rect_coords)], 
         ybottom = min(data_format[indices_i,c("hr", "map")]), 
         xright = rect_coords$t[-1], 
         ytop = max(data_format[indices_i,c("hr", "map")]),
         col = col_vec[-nrow(rect_coords)],
         border = NA)
    points(t_grid, data_format[indices_i,"hr"], col = 'aquamarine')
    points(t_grid, data_format[indices_i,"map"], col = 'orange')
    
    # Vitals 1 & 4 -------------------------------------------------------------
    plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=c(min(data_format[indices_i,c("hemo", "lactate")]),
                                                    max(data_format[indices_i,c("hemo", "lactate")])), 
         main=paste0("EID = ", i, " truth"), xlab='time', ylab=NA, xaxt='n', 
         col.main='green', col.axis='green')
    
    rect(xleft = rect_coords$t[-nrow(rect_coords)], 
         ybottom = min(data_format[indices_i,c("hemo", "lactate")]), 
         xright = rect_coords$t[-1], 
         ytop = max(data_format[indices_i,c("hemo", "lactate")]),
         col = col_vec[-nrow(rect_coords)],
         border = NA)
    points(t_grid, data_format[indices_i,"hemo"], col = 'white')
    points(t_grid, data_format[indices_i,"lactate"], col = 'purple')
    
    # BAR PLOTS --------------------------------------------------------------
    barplot(state_results[[seed_focus]][1:S, indices_i],
            col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
            xlab='time', space=0, col.main='green', border=NA,
            xlim=range(pb) + c(-0.5,0.5), main = paste0("EID = ", i, " fit"))
    grid( nx=20, NULL, col='white')
    legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'Baseline', 'State 2', 'State 3'),
            pch=15, pt.cex=1.5,
            col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'))
    axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
    axis( side=2, at=0:1, col.axis='green')
}
dev.off()
