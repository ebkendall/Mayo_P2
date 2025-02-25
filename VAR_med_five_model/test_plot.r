library(matrixStats)
library(plotrix)

# Function to change transparency of colors # ----------------------------------
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

B_curr = do.call('c', B)
Hr_chain = Y[,"hr"]
Map_chain = Y[,"map"]
Hc_chain = Y[,"hemo"]
La_chain = Y[,"lactate"]

pdf(paste0('Plots/initialized_Y.pdf'))
panel_dim = c(4,1)
inset_dim = c(0,-.18)
par(mfrow=panel_dim, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
    if(which(EIDs == i) %% 100 == 0) print(which(EIDs == i))
    
    indices_i = (data_format[,'EID']==i)
    n_i = sum(indices_i)
    
    t_grid = round(data_format[indices_i, 'time'] / 60, digits = 3)
    t_grid_bar = 1:length(t_grid)
    rbc_times_bar = which(data_format[data_format[,'EID']==i, 'RBC_ordered'] != 0)
    if(simulation) {
        rbc_admin = c(head(data_format[data_format[,'EID']==i, "n_RBC_admin"], 1),
                      diff(data_format[data_format[,'EID']==i, "n_RBC_admin"]))
        rbc_admin_times_bar = which(rbc_admin != 0)
    } else {
        rbc_admin_times_bar = which(data_format[data_format[,'EID']==i, 'RBC_admin'] != 0)   
    }
    rbc_times = t_grid[rbc_times_bar]
    rbc_admin_times = t_grid[rbc_admin_times_bar]
    
    # Put this on the correct scale as the t_grid
    b_i = B_curr[indices_i]
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
    
    pb = barplot(matrix(1/5, nrow = 5, ncol = n_i), 
                 col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                 xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
    
    # Heart Rate and MAP double plot -----------------------------------------
    if(mean(data_format[indices_i, 'clinic_rule']) != 0) {
        title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                            mean(data_format[indices_i, 'RBC_rule']),
                            ', clinic = ', mean(data_format[indices_i, 'clinic_rule']))
    } else {
        title_name = paste0('Heart Rate & MAP: ', i, ', RBC Rule = ', 
                            mean(data_format[indices_i, 'RBC_rule']))
    }
    
    hr_map_ylim = c(min(Hr_chain[indices_i], Map_chain[indices_i]), 
                    max(Hr_chain[indices_i], Map_chain[indices_i]))
    
    # Make a new plot to add the background color
    plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
         xlab='time', ylab=NA, xaxt='n', col.main='green',
         col.axis='green')
    
    rect(xleft = rect_coords$t[-nrow(rect_coords)], 
         ybottom = hr_map_ylim[1], 
         xright = rect_coords$t[-1], 
         ytop = hr_map_ylim[2],
         col = col_vec[-nrow(rect_coords)],
         border = NA)
    
    points( x = 1:n_i, y=Hr_chain[indices_i], pch=c(20 + otype[indices_i, 2]), col = 'aquamarine') 
    points( x = 1:n_i, y=Map_chain[indices_i], pch=c(20 + otype[indices_i, 3]), col = 'orange') 
    legend( 'topright', inset=c(0,-.18), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'HR', 'MAP'), pch=15, pt.cex=1.5, 
            col=c( 'aquamarine', 'orange'))
    
    
    # Hemoglobin and Lactate double plot -------------------------------------
    if(mean(data_format[indices_i, 'clinic_rule']) != 0) {
        title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ', 
                            mean(data_format[indices_i, 'RBC_rule']),
                            ', clinic = ', mean(data_format[indices_i, 'clinic_rule']))
    } else {
        title_name = paste0('Hemoglobin & Lactate: ', i, ', RBC Rule = ',
                            mean(data_format[indices_i, 'RBC_rule']))
    }
    
    hr_map_ylim = c(min(Hc_chain[indices_i], La_chain[indices_i]), 
                    max(Hc_chain[indices_i], La_chain[indices_i]))
    
    plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
         xlab='time', ylab=NA, xaxt='n', col.main='green',
         col.axis='green')
    
    rect(xleft = rect_coords$t[-nrow(rect_coords)], 
         ybottom = hr_map_ylim[1], 
         xright = rect_coords$t[-1], 
         ytop = hr_map_ylim[2],
         col = col_vec[-nrow(rect_coords)],
         border = NA)
    
    points( x = 1:n_i, y=Hc_chain[indices_i], pch=c(20 + otype[indices_i, 1]), col = 'aquamarine') 
    points( x = 1:n_i, y=La_chain[indices_i], pch=c(20 + otype[indices_i, 4]), col = 'orange') 
    legend( 'topright', inset=c(0,-.18), xpd=T, horiz=T, bty='n', x.intersp=.75,
            legend=c( 'HEMO', 'LACT'), pch=15, pt.cex=1.5, 
            col=c( 'aquamarine', 'orange'))
}
dev.off()
