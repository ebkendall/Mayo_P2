load('Data/long_data_agg.rda') # long_data_agg
load('Data/all_keys.rda')      # all_keys
load('Data/cov_info.rda')      # cov_info
load('Data/timing_issues.rda') # timing_issues

times = rep(NA, length(long_data_agg))
for(i in 1:length(long_data_agg)) {
    times[i] = nrow(long_data_agg[[i]]$covariates) / 4
}

clinic_review = read.csv('../selected_encounters_update.csv')
clinical_id = as.numeric(clinic_review$key)
clinical_id_yes = clinical_id[clinic_review$Bleeding.event. %in% c("Yes", "Y")]
clinical_id_no = clinical_id[clinic_review$Bleeding.event. %in% c("No", "N")]

# Temporarily sub-setting long_data_agg ----------------------------------------
long_data_agg_sub = list()
all_keys_temp = NULL
ldas = 1
for(i in 1:length(long_data_agg)) {
    # Removing time issues
    if(long_data_agg[[i]]$time_flag == 0) {
        # Need at least 2 observations
        if(nrow(long_data_agg[[i]]$covariates) > 1) {
            
            long_data_agg_sub[[ldas]] = long_data_agg[[i]]
            all_keys_temp = c(all_keys_temp, long_data_agg[[i]]$key)
            
            if(!(long_data_agg[[i]]$key %in% long_data_agg_sub[[ldas]]$covariates[,'EID'])) {print("wrong ID")}
            
            ldas = ldas + 1    
        }
    } else {
        print(paste0("time issue: ", long_data_agg[[i]]$key))
    }
}
for(i in 1:length(long_data_agg_sub)) {
    if(long_data_agg_sub[[i]]$key != all_keys_temp[i]) print("wrong order")
}
rm(long_data_agg)
rm(cov_info)

# Check the level-of-care for each patient -------------------------------------
level_of_care = read.csv('raw_data/jw_patient_level_of_care.csv')
care_props = rep(0, length(long_data_agg_sub))
total_time_icu = rep(0, length(long_data_agg_sub))
min_max_icu_time = matrix(-1, nrow = length(long_data_agg_sub), ncol = 2)

print("Going through level of care to find ICU times")
for(i in 1:length(long_data_agg_sub)) {
    ii = all_keys_temp[i]
    temp = long_data_agg_sub[[i]]$covariates
    if(ii != long_data_agg_sub[[i]]$key[1]) print(paste0(i, ": wrong ID!"))
    
    care_time_df = matrix(nrow = nrow(temp), ncol = 3)
    care_time_df[,1] = temp[,'EID']
    care_time_df[,2] = temp[,'time']
    
    level_of_care_patient = level_of_care[level_of_care$key %in% temp[,"EID"],,drop=F]
    level_of_care_patient$level_of_care_datetime = level_of_care_patient$level_of_care_datetime / 60
    
    sub_level = level_of_care_patient[level_of_care_patient$key == ii, ,drop=F]
    sub_df    = temp[temp[,"EID"] == ii, ,drop=F]
    sub_care  = care_time_df[care_time_df[,1] == ii, ,drop=F]
    level_times = sub_level$level_of_care_datetime
    
    if(sum(sub_level$patient_level_of_care == "Intensive Care") != nrow(sub_level)) {
        for(j in 1:nrow(sub_df)) {
            time_j = sub_df[j,"time"]
            
            if(time_j < level_times[1]) {
                ind = which(level_times == level_times[1])
                sub_care[j,3] = sub_level$patient_level_of_care[ind]
            } else if(time_j >= tail(level_times,1)) {
                ind = which(level_times == tail(level_times,1))
                sub_care[j,3] = sub_level$patient_level_of_care[ind]
            } else {
                end_pts = c(max(level_times[level_times <= time_j]), 
                            min(level_times[level_times >= time_j]))
                
                if(sum(is.finite(end_pts)) != length(end_pts)) print(paste0("issue: ", i))
                
                ind = which(level_times == end_pts[1])
                care_ind = NULL
                if("Intensive Care" %in% sub_level$patient_level_of_care[ind] &
                   sub_level$patient_level_of_care[max(ind) + 1] == "Intensive Care") {
                    care_ind = ind[which("Intensive Care" %in% sub_level$patient_level_of_care[ind])]
                } else {
                    care_ind = max(ind)
                }
                sub_care[j,3] = sub_level$patient_level_of_care[care_ind]
            }
        }
        care_time_df[, 3] = sub_care[,3]   
    } else {
        care_time_df[, 3] = rep("Intensive Care", nrow(care_time_df))
        sub_care[,3] = rep("Intensive Care", nrow(sub_care))
    }
    
    care_props[i] = sum(sub_care[,3] == "Intensive Care") / nrow(sub_care)
    
    long_data_agg_sub[[i]]$covariates = cbind(long_data_agg_sub[[i]]$covariates, 
                                              sub_care[,3])
    
    if(care_props[i] > 0) {
        min_icu = min(which(sub_care[,3] == "Intensive Care"))
        max_icu = max(which(sub_care[,3] == "Intensive Care"))
        total_time_icu[i] = sum(sub_care[,3] == "Intensive Care")
        min_max_icu_time[i,] = c(min_icu, max_icu)
    }
}
print("Done")

print(paste0("Remaining clinical patients: ", sum(clinical_id %in% all_keys_temp),
             " of the initial ", length(clinical_id)))

# Making long_data_agg cleaner -------------------------------------------------
long_data_agg_sub2 = vector(mode = 'list', length = length(all_keys_temp))
ldas = 1
for(i in 1:length(long_data_agg_sub2)) {
    long_data_agg_sub2[[ldas]] = long_data_agg_sub[[i]]$covariates
    ldas = ldas + 1
}
long_data_agg_sub = long_data_agg_sub2
rm(long_data_agg_sub2)

# Formatting the existing data into one data set -------------------------------
print("Compiling all information into data_format")
print(paste0("Number of subjects: ", length(all_keys_temp)))
data_format = do.call('rbind', long_data_agg_sub)

data_format = cbind(data_format, rep(0, nrow(data_format)), rep(0, nrow(data_format)))
colnames(data_format) = c('EID', 'time', 'temperature', 'hemo', 'map', 'hr', 'lactate', 'RBC',
                          'n_labs', 'n_RBC', 'levelCare', 'RBC_rule', 'clinic_rule')
# ------------------------------------------------------------------------------
# (2) Add the new RBC transfusions ---------------------------------------------
# ------------------------------------------------------------------------------
transfusions = read.csv('raw_data/jw_transfusions.csv')
transfusions = transfusions[transfusions$OrderedProduct == "RBC", ]

times = transfusions[,c("key",
                        "transfus_order_datetime", 
                        "transfus_admin_start_datetime", 
                        "OrderedProduct", 
                        "Order_Description")]

data_keys = unique(data_format[,"EID"])
times = times[times$key %in% data_keys, ]
times = times[times$transfus_order_datetime > 0, ]
times$transfus_order_datetime = times$transfus_order_datetime / 60
times$transfus_admin_start_datetime = times$transfus_admin_start_datetime / 60
times = times[!is.na(times$key), ]

RBC_ordered = 0
RBC_admin = 0
data_format = cbind(cbind(data_format, RBC_ordered), RBC_admin)

level_of_care_vec = data_format[,"levelCare"]
data_format = data_format[,(colnames(data_format) != 'levelCare')]
data_format = apply(data_format, 2, as.numeric)

print("Adding transfusion order and admin time")
for(i in 1:length(data_keys)) {
    if(data_keys[i] %in% times$key) {
        sub_dat = data_format[data_format[,"EID"] == as.numeric(data_keys[i]), ,drop = F]
        sub_rbc = times[times$key == data_keys[i], , drop = F]
        for(j in 1:(nrow(sub_dat) - 1)) {
            timed_rbc_order = which(sub_rbc$transfus_order_datetime > sub_dat[j,"time"] & sub_rbc$transfus_order_datetime <= sub_dat[j+1,"time"])
            if(length(timed_rbc_order) > 0) {
                sub_dat[j+1, "RBC_ordered"] = length(timed_rbc_order)
            }
            
            timed_rbc_admin = which(sub_rbc$transfus_admin_start_datetime > sub_dat[j,"time"] & sub_rbc$transfus_admin_start_datetime <= sub_dat[j+1,"time"])
            if(length(timed_rbc_admin) > 0) {
                sub_dat[j+1, "RBC_admin"] = length(timed_rbc_admin)
            }
        }
        
        data_format[data_format[,"EID"] == as.numeric(data_keys[i]), ] = sub_dat
    }
}
print("Done")

# Adding n_RBC_ordered and n_RBC_admin
n_RBC_ordered = 0
n_RBC_admin = 0
data_format = cbind(cbind(data_format, n_RBC_ordered), n_RBC_admin)

print("Calculating n_RBC_ordered and n_RBC_admin")
for(i in unique(data_format[,'EID'])){
    data_format[data_format[,'EID'] == i, 'n_RBC_ordered'] = cumsum(data_format[data_format[,'EID'] == i, 'RBC_ordered'])
    data_format[data_format[,'EID'] == i, 'n_RBC_admin'] = cumsum(data_format[data_format[,'EID'] == i, 'RBC_admin'])
}
print("Done")

# Adding the clinical rule -----------------------------------------------------
data_format[data_format[,"EID"] %in% clinical_id_yes, "clinic_rule"] = 1 # bleed
data_format[data_format[,"EID"] %in% clinical_id_no, "clinic_rule"] = -1 # no bleed

# Remove any erroneous data points ---------------------------------------------
for(i in 1:nrow(data_format)) {
    # hemo
    if(!is.na(data_format[i,'hemo'])) {
        if(data_format[i,'hemo'] < 0 | data_format[i,'hemo'] > 20) {
            data_format[i, 'hemo'] = NA
        } 
    }
    # hr
    if(!is.na(data_format[i,'hr'])) {
        if(data_format[i,'hr'] < 20 | data_format[i,'hr'] > 200) {
            data_format[i, 'hr'] = NA
        } 
    }
    # map
    if(!is.na(data_format[i,'map'])) {
        if(data_format[i,'map'] < 25 | data_format[i,'map'] > 150) {
            data_format[i, 'map'] = NA
        } 
    }
    
}

print(paste0("Remaining clinical patients: ", sum(clinical_id %in% data_format[,"EID"]),
             " of the initial ", length(clinical_id)))

if(length(level_of_care_vec) != nrow(data_format)) {print("level of care vec wrong")}

save(data_format, file = "Data/data_format.rda")
save(level_of_care_vec, file = "Data/level_of_care.rda")