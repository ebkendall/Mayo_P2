args = commandArgs(TRUE)
sampling_num = as.numeric(args[1])
index_seeds = 1:3
it_num = 2

true_par = c(2, -2, 0.4054651, 0.4054651)
par_index = list()
par_index$mu = 1:2
par_index$t_p = 3:4

labels = c("mu1", "mu2", "logit baseline 1 -> 2", "logit baseline 2 -> 1") 

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = length(index_seeds))

ind = 0

for(seed in index_seeds){
    
    it_seq = 1:it_num
    ind = ind + 1
    
    for(it in it_seq) {
        file_name = paste0('Model_out/mcmc_out_',toString(seed),'_', 'it', 
                           it, '_samp', sampling_num, '_sim.rda') 
        
        load(file_name)
        print(paste0(ind, ": ", file_name))
        print("accept")
        print(mcmc_out$accept)
        print("pscale")
        print(mcmc_out$pscale)
        
        par_index = mcmc_out$par_index
        
        if(it == 1) {
            chain_list[[ind]] = mcmc_out$chain[1000:2000,]
        } else {
            chain_list[[ind]] = rbind(chain_list[[ind]], mcmc_out$chain)
        }
        
        rm(mcmc_out)
    }
}

stacked_chains = do.call( rbind, chain_list)

pdf_title = paste0('trace_plot_samp', sampling_num, '.pdf')
pdf(pdf_title)
par(mfrow=c(3, 2))
lab_ind = 0
for(s in names(par_index)){
    temp_par = par_index[[s]]
    for(r in temp_par){
        # lab_ind = lab_ind + 1
        lab_ind = r
        parMean = round( mean(stacked_chains[,r]), 4)
        parMedian = round( median(stacked_chains[,r]), 4)
        upper = quantile( stacked_chains[,r], prob=.975)
        lower = quantile( stacked_chains[,r], prob=.025)
        
        title_color = "black"
        
        y_limit = range(stacked_chains[,r])
        plot( NULL, ylab=NA, main=labels[lab_ind], xlim=c(1,nrow(chain_list[[1]])),
              ylim=y_limit, xlab = paste0("95% CI: [", round(lower, 4),
                                          ", ", round(upper, 4), "]"),
              col.main = title_color)
        
        for(seed in 1:length(chain_list)) { lines( chain_list[[seed]][,r], type='l', col=seed) }
        
        x_label = paste0('Mean =',toString(parMean),
                         ' Median =',toString(parMedian),
                         ' True =', round(true_par[r], 3))
        
        hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA, freq=FALSE,
              xlab=x_label)
        abline( v=upper, col='red', lwd=2, lty=2)
        abline( v=lower, col='purple', lwd=2, lty=2)
        abline( v=true_par[r], col='green', lwd=2, lty=2)
    }   
}


dev.off()
