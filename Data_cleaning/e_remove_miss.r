# Remove long strands of missingness in the real data

load_names_df = c('Data/data_format_train.rda', 'Data/data_format_train_large.rda')
load_names_Dn = c('Data/Dn_omega.rda', 'Data/Dn_omega_large.rda')

save_names_df = c('Data/data_format_train_miss.rda', 'Data/data_format_train_miss_large.rda')
save_names_Dn = c('Data/Dn_omega_miss.rda', 'Data/Dn_omega_miss_large.rda')

for(df_name in 1:2) {
    
    load(load_names_df[df_name])
    load(load_names_Dn[df_name])
    
    EIDs = unique(data_format[,"EID"])
    
    first_and_last_ind = matrix(nrow = length(EIDs), ncol = 2)
    
    for(ii in 1:length(EIDs)) {
        
        sub_dat = data_format[data_format[,"EID"] == EIDs[ii], ]
        
        first_ind = 1
        last_ind = nrow(sub_dat)
        
        # Check for missingness at *start* of encounter
        if(is.na(sub_dat[1, "hr"]) || is.na(sub_dat[1,"map"])) {
            missing_hr_sum = cumsum(is.na(sub_dat[, "hr"]))
            missing_map_sum = cumsum(is.na(sub_dat[, "map"]))
            
            indexing = 1:nrow(sub_dat)
            
            diff_hr = abs(indexing - missing_hr_sum)
            diff_map = abs(indexing - missing_map_sum)
            
            first_hr_obs = 1
            first_map_obs = 1
            if(0 %in% diff_hr) { first_hr_obs = max(which(diff_hr == 0)) + 1 }
            if(0 %in% diff_map){ first_map_obs = max(which(diff_map == 0)) + 1 }
            
            first_ind = max(c(first_hr_obs, first_map_obs))
        }
        
        # Check for missingness at *end* of encounter
        if(is.na(sub_dat[nrow(sub_dat), "hr"]) || is.na(sub_dat[nrow(sub_dat),"map"])) {
            
            # Flip sub_dat "upside down"
            sub_dat_flip = sub_dat[nrow(sub_dat):1, ]
            
            missing_hr_sum = cumsum(is.na(sub_dat_flip[, "hr"]))
            missing_map_sum = cumsum(is.na(sub_dat_flip[, "map"]))
            
            indexing = 1:nrow(sub_dat)
            
            diff_hr = abs(indexing - missing_hr_sum)
            diff_map = abs(indexing - missing_map_sum)
            
            first_hr_obs = 1
            first_map_obs = 1
            if(0 %in% diff_hr) { first_hr_obs = max(which(diff_hr == 0)) + 1 }
            if(0 %in% diff_map){ first_map_obs = max(which(diff_map == 0)) + 1 }
            
            missing_last_ind = max(c(first_hr_obs, first_map_obs))
            
            last_ind = nrow(sub_dat) - missing_last_ind + 1
        }
        
        first_and_last_ind[ii, ] = c(first_ind, last_ind)
    }
    
    new_data_format = NULL
    new_Dn_omega = list()
    
    for(ii in 1:length(EIDs)) {
        
        print(ii)
        
        index_keep = first_and_last_ind[ii,1]:first_and_last_ind[ii,2]
        sub_dat = data_format[data_format[,"EID"] == EIDs[ii], ]
        
        new_data_format = rbind(new_data_format, sub_dat[index_keep, ])
        
        new_Dn_omega[[ii]] = Dn_omega[[ii]][index_keep]
    }
    
    data_format = new_data_format
    Dn_omega = new_Dn_omega
    
    save(data_format, file = save_names_df[df_name])
    save(Dn_omega, file = save_names_Dn[df_name]) 
    
    rm(data_format)
    rm(Dn_omega)
    rm(new_data_format)
    rm(new_Dn_omega)
}
