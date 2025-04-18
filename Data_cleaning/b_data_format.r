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

clinical_lab = clinical_id
# clinical_lab = c(109125, 111375, 133750, 156325, 165725, 195475, 198975, 208100, 327375, 360525,
#                  431400, 467200, 494300, 531650, 533825, 543625, 588100, 622450, 633100, 697600,
#                  727675, 750900, 758025, 781875, 785950, 801300, 820775, 827350, 841925, 843000)

# ------------------------------------------------------------------------------
# (1) Filtering people based on undesirable characteristics --------------------
# ------------------------------------------------------------------------------

# Removing time issues
all_keys_temp = all_keys[!(all_keys %in% timing_issues)]

# Removing patients who died in the ICU
deaths = cov_info[cov_info[,"icu_death"] == 1, "key"]
all_keys_temp = all_keys_temp[!(all_keys_temp %in% deaths)]

# Isolating focus to patients who have had a long enough stay ------------------
max_times = matrix(data = 0, nrow = length(all_keys_temp), ncol = 2)
max_times[,1] = all_keys_temp

for(i in 1:length(long_data_agg)) {
    if(long_data_agg[[i]]$key %in% all_keys_temp) {
        if(nrow(long_data_agg[[i]]$covariates) > 1) {
            max_times[max_times[,1] == long_data_agg[[i]]$key, 2] = 
                nrow(long_data_agg[[i]]$covariates)
        }
    }
}

# Only take patients with >= 6hr length of stay (>= 24 rows)
all_keys_temp = max_times[max_times[,2] >= 24,1]

# Temporarily sub-setting long_data_agg ----------------------------------------
long_data_agg_sub = vector(mode = 'list', length = length(all_keys_temp))
ldas = 1
for(i in 1:length(long_data_agg)) {
    if(long_data_agg[[i]]$key %in% all_keys_temp) {
        long_data_agg_sub[[ldas]] = long_data_agg[[i]]
        
        # Only take the first 48 hours (192 rows) of a patient encounter -------
        if(nrow(long_data_agg_sub[[ldas]]$covariates) > 192) {
            long_data_agg_sub[[ldas]]$covariates = long_data_agg_sub[[ldas]]$covariates[1:192, ,drop=F]
        }
        ldas = ldas + 1
    }
}
for(i in 1:length(long_data_agg_sub)) {
    if(long_data_agg_sub[[i]]$key != all_keys_temp[i]) print("wrong order")
}
rm(long_data_agg)
rm(cov_info)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Go through and look at the missingness of the data and see if it is because of 
# level of care
level_of_care = read.csv('raw_data/jw_patient_level_of_care.csv')
care_props = rep(0, length(long_data_agg_sub))
total_time_icu = rep(0, length(long_data_agg_sub))
min_max_icu_time = matrix(-1, nrow = length(long_data_agg_sub), ncol = 2)

check_patients_ind = c(5966, 7321, 12557, 21995)

print("Going through level of care to find ICU times")
for(i in 1:length(long_data_agg_sub)) {
    ii = all_keys_temp[i]
    temp = long_data_agg_sub[[i]]$covariates
    if(ii != long_data_agg_sub[[i]]$key[1]) print(paste0(i, ": wrong ID!"))
    
    care_time_df = matrix(nrow = nrow(temp), ncol = 3)
    care_time_df[,1] = temp[,'EID']
    care_time_df[,2] = temp[,'time']
    
    level_of_care_patient = level_of_care[level_of_care$key %in% temp[,"EID"],]
    level_of_care_patient$level_of_care_datetime = level_of_care_patient$level_of_care_datetime / 60
    
    sub_level = level_of_care_patient[level_of_care_patient$key == ii, ]
    sub_df    = temp[temp[,"EID"] == ii, ]
    sub_care  = care_time_df[care_time_df[,1] == ii, ]
    level_times = sub_level$level_of_care_datetime
    
    if(sum(sub_level$patient_level_of_care == "Intensive Care") != nrow(sub_level)) {
        for(j in 1:nrow(sub_df)) {
            time_j = sub_df[j,"time"]
            
            if(time_j < level_times[1]) {
                ind = which(level_times == level_times[1])
                sub_care[j,3] = sub_level$patient_level_of_care[ind]
            } else if(time_j > tail(level_times,1)) {
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
        # long_data_agg_sub[[i]]$covariates = long_data_agg_sub[[i]]$covariates[min_icu:max_icu, ,drop=F]   
    }
}
print("Done")

# Remove subjects with 0 time in intensive level of care -----------------------
all_keys_temp2 = all_keys_temp[total_time_icu > 0]
long_data_agg_sub2 = vector(mode = 'list', length = length(all_keys_temp2))
ldas = 1
for(i in 1:length(long_data_agg_sub)) {
    if(long_data_agg_sub[[i]]$key %in% all_keys_temp2) {
        long_data_agg_sub2[[ldas]] = long_data_agg_sub[[i]]
        ldas = ldas + 1
    }
}
for(i in 1:length(long_data_agg_sub2)) {
    if(long_data_agg_sub2[[i]]$key != all_keys_temp2[i]) print("wrong order")
}

all_keys_temp = all_keys_temp2
long_data_agg_sub = long_data_agg_sub2
rm(all_keys_temp2)
rm(long_data_agg_sub2)

# Choose subjects with >= 50% HR and MAP measurements --------------------------
enough_dat = rep(0, length(all_keys_temp))
for(i in 1:length(enough_dat)) {
    temp = long_data_agg_sub[[i]]$covariates
    hr_ratio = sum(!is.na(temp[,"hr"])) / nrow(temp)
    map_ratio = sum(!is.na(temp[,"map"])) / nrow(temp)
    
    if(hr_ratio > 0.5 & map_ratio > 0.5) {
        enough_dat[i] = 1
    }
}

all_keys_temp = all_keys_temp[enough_dat == 1]

print(paste0("Remaining clinical patients: ", sum(clinical_lab %in% all_keys_temp),
             " of the initial ", length(clinical_id)))

# # Adding back all clinically annotated patients --------------------------------
# all_keys_temp = unique(c(all_keys_temp, clinical_id))
# print(paste0("Adding back clinical patients: ", sum(clinical_lab %in% all_keys_temp),
#              " of the initial ", length(clinical_id)))

# Formatting the existing data into one data set -------------------------------
print("Compiling all information into data_format")
print(paste0("Number of subjects: ", length(all_keys_temp)))

data_format = NULL
for(i in 1:length(long_data_agg_sub)) {
    if(long_data_agg_sub[[i]]$key %in% all_keys_temp) {
        data_format = rbind(data_format, long_data_agg_sub[[i]]$covariates)
    }
}
print("Done")

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

print("Adding transfusion order and admin time")
for(i in 1:length(data_keys)) {
    if(data_keys[i] %in% times$key) {
        sub_dat = data_format[data_format[,"EID"] == data_keys[i], ,drop = F]
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
        
        data_format[data_format[,"EID"] == data_keys[i], ] = sub_dat
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
# ------------------------------------------------------------------------------
# (2) Adding new variables such as RBC rule ------------------------------------
# ------------------------------------------------------------------------------

level_of_care_vec = data_format[,11]
data_format = data_format[,-11]
data_format = apply(data_format, 2, as.numeric)

# Get patients that had at least 3 RBC at some point in their stay
min_three_RBC = unique(data_format[data_format[,"n_RBC_admin"] >= 3, "EID"])

bleed_pat = min_three_RBC

# Adding the rule that a patient is bleeding if >= 3 in 12hrs or >= 6 in 24hrs
for(i in 1:length(bleed_pat)) {
    sub_dat = data_format[data_format[,"EID"] == bleed_pat[i], ]
    
    # Check in any 12 hour period
    max_time = tail(sub_dat[,"time"], 1)
    when_rbc = c(1, which(diff(sub_dat[,"n_RBC_admin"]) != 0))
    
    for(j in 1:length(when_rbc)) {
        s_time = sub_dat[when_rbc[j], "time"]
        e_time_12 = s_time + 720
        e_time_24 = s_time + 1440
        RBC_diff_12 = RBC_diff_24 = 0
        
        if (e_time_12 <= max_time) {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            ind_12 = order(abs(sub_dat[,"time"] - e_time_12))[1]
            RBC_diff_12 = sub_dat[ind_12, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        } else {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
            RBC_diff_12 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        }
        if (e_time_24 <= max_time) {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            ind_24 = order(abs(sub_dat[,"time"] - e_time_24))[1]
            RBC_diff_24 = sub_dat[ind_24, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        } else {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
            RBC_diff_24 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        }
        
        if(RBC_diff_12 >=3 | RBC_diff_24 >= 6) {
            data_format[data_format[,"EID"] == bleed_pat[i], "RBC_rule"] = 1
            break
        }
        
    }
    
}

# Adding the clinical rule
data_format[data_format[,"EID"] %in% clinical_id_yes, "clinic_rule"] = 1 # bleed
data_format[data_format[,"EID"] %in% clinical_id_no, "clinic_rule"] = -1 # no bleed

# ------------------------------------------------------------------------------
# Removing pacing patients based on own heuristic ------------------------------
# ------------------------------------------------------------------------------
pace_ids = NULL
for(i in unique(data_format[,"EID"])) {
    sub_id_hr = data_format[data_format[,"EID"] == i, "hr"]
    diff_hr = diff(sub_id_hr)
    conseq_same = which(diff_hr == 0)
    if(sum(diff(conseq_same) == 1) > 6) {
        pace_ids = c(pace_ids, i)
    }
}

level_of_care_vec = level_of_care_vec[!(data_format[,"EID"] %in% pace_ids)]
data_format = data_format[!(data_format[,"EID"] %in% pace_ids), ]

# ------------------------------------------------------------------------------
# (5) Remove any eroneous data points ------------------------------------------
# ------------------------------------------------------------------------------
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

print(paste0("Remaining clinical patients: ", sum(clinical_lab %in% data_format[,"EID"]),
             " of the initial ", length(clinical_id)))

if(length(level_of_care_vec) != nrow(data_format)) {print("level of care vec wrong")}

save(data_format, file = "Data/data_format.rda")
save(level_of_care_vec, file = "Data/level_of_care.rda")



# Want subjects with at least 6 hours of intensive care attention --------------
# # all_keys_temp2 = all_keys_temp[care_props != 0]
# all_keys_temp2 = all_keys_temp[total_time_icu >= 24]
# 
# long_data_agg_sub2 = vector(mode = 'list', length = length(all_keys_temp2))
# ldas = 1
# for(i in 1:length(long_data_agg_sub)) {
#     if(long_data_agg_sub[[i]]$key %in% all_keys_temp2) {
#         long_data_agg_sub2[[ldas]] = long_data_agg_sub[[i]]
#         ldas = ldas + 1
#     }
# }
# for(i in 1:length(long_data_agg_sub2)) {
#     if(long_data_agg_sub2[[i]]$key != all_keys_temp2[i]) print("wrong order")
# }
# 
# all_keys_temp = all_keys_temp2
# long_data_agg_sub = long_data_agg_sub2
# 
# print(paste0("Remaining clinical patients: ", sum(clinical_lab %in% all_keys_temp),
#              " of the initial ", length(clinical_id)))
# 
# rm(all_keys_temp2)
# rm(long_data_agg_sub2)
# # -----------------------------------------------------------------------------

# # Restricting the amount of information we have on each patient to 
# print("Counting how many observations each patient has")
# unique_id = unique(data_format[,"EID"])
# total_obs = table(as.numeric(data_format[,"EID"]))
# total_obs_names = names(total_obs)
# 
# if(sum(unique_id != total_obs_names) != 0) {print("mismatch with table")}
# 
# id_keep = unique_id[total_obs >= 40]
# print(paste0("Remaining clinical patients: ", sum(clinical_lab %in% id_keep)))
# data_format = data_format[data_format[,"EID"] %in% id_keep, ]

