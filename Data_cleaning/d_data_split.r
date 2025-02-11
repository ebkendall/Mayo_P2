# Load in the large format data ------------------------------------------------
load('Data/data_format.rda')
load('Data/level_of_care.rda')
load('Data/cov_info.rda')  

EIDs = unique(data_format[,"EID"])

# ------------------------------------------------------------------------------
# FINAL FORMATTING TO DETERMINE TRAINING SET -----------------------------------
# ------------------------------------------------------------------------------

# Removing patients who died in the ICU  ---------------------------------------
alive = cov_info[cov_info[,"icu_death"] != 1, "key"]
ind_alive = which(data_format[,"EID"] %in% alive)
data_format = data_format[ind_alive, ]
level_of_care_vec = level_of_care_vec[ind_alive]
EIDs_sub = unique(data_format[,"EID"])

# Only take patients with >= 6hr length of stay (>= 24 rows) -------------------
patient_stay_time = table(data_format[,"EID"])
if(sum(names(patient_stay_time) == EIDs_sub) != length(EIDs_sub)) print("error")
EIDs_sub = EIDs_sub[patient_stay_time >= 24]
ind_min_6 = which(data_format[,"EID"] %in% EIDs_sub)
data_format = data_format[ind_min_6, ]
level_of_care_vec = level_of_care_vec[ind_min_6]

if(nrow(data_format) != length(level_of_care_vec)) print("error")

# Only take the first 48 hours (192 rows) of a patient encounter ---------------
patient_stay_time = table(data_format[,"EID"])
if(sum(names(patient_stay_time) == EIDs_sub) != length(EIDs_sub)) print("error")
EID_more_48 = EIDs_sub[patient_stay_time > 192]

for(i in 1:length(EID_more_48)) {
    print(i)
    id_i = EID_more_48[i]
    ind_i = which(data_format[,"EID"] == id_i)
    remove_inds = ind_i[193:length(ind_i)]
    data_format = data_format[-remove_inds, ]
    level_of_care_vec = level_of_care_vec[-remove_inds]
}

if(nrow(data_format) != length(level_of_care_vec)) print("error")
EIDs_sub = unique(data_format[,"EID"])

# Remove subjects with 0 time in intensive level of care -----------------------
EID_no_ICU = NULL
for(i in EIDs_sub) {
    sub_df = data_format[data_format[,"EID"] == i, ]
    sub_loc = level_of_care_vec[data_format[,"EID"] == i]
    if(sum(sub_loc == "Intensive Care") == 0) {
        print(i)
        EID_no_ICU = c(EID_no_ICU, i)
    }
}

ind_ICU = which(!(data_format[,"EID"] %in% EID_no_ICU))
data_format = data_format[ind_ICU, ]
level_of_care_vec = level_of_care_vec[ind_ICU]
if(nrow(data_format) != length(level_of_care_vec)) print("error")
EIDs_sub = unique(data_format[,"EID"])

# Choose subjects with >= 50% HR and MAP measurements --------------------------
enough_dat = rep(0, length(EIDs_sub))
for(i in 1:length(EIDs_sub)) {
    if(i %% 100 == 0) print(i)
    sub_df = data_format[data_format[,"EID"] == EIDs_sub[i], ]
    
    hr_ratio = sum(!is.na(sub_df[,"hr"])) / nrow(sub_df)
    map_ratio = sum(!is.na(sub_df[,"map"])) / nrow(sub_df)

    if(hr_ratio >= 0.5 & map_ratio >= 0.5) {
        enough_dat[i] = 1
    }
}

EIDs_sub = EIDs_sub[enough_dat == 1]

ind_dat = which(data_format[,"EID"] %in% EIDs_sub)
data_format = data_format[ind_dat, ]
level_of_care_vec = level_of_care_vec[ind_dat]
if(nrow(data_format) != length(level_of_care_vec)) print("error")

# Removing pacing patients based on own heuristic ------------------------------
pace_ids = NULL
for(i in EIDs_sub) {
    if(i %% 100 == 0) print(i)
    sub_id_hr = data_format[data_format[,"EID"] == i, "hr"]
    diff_hr = diff(sub_id_hr)
    conseq_same = which(diff_hr == 0)
    if(sum(diff(conseq_same) == 1) > 6) {
        pace_ids = c(pace_ids, i)
    }
}

ind_no_pace = which(!(data_format[,"EID"] %in% pace_ids))
data_format = data_format[ind_no_pace, ]
level_of_care_vec = level_of_care_vec[ind_no_pace]
if(nrow(data_format) != length(level_of_care_vec)) print("error")
EIDs_sub = unique(data_format[,"EID"])

# Adding RBC rule --------------------------------------------------------------

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

# Remove patients who's medication records had inconsistencies -----------------
load('Data/flagged_patients_med.rda')
ind_no_flag = which(!(data_format[,"EID"] %in% flagged_patients))
data_format = data_format[ind_no_flag, ]
level_of_care_vec = level_of_care_vec[ind_no_flag]
if(nrow(data_format) != length(level_of_care_vec)) print("error")
EIDs_sub = unique(data_format[,"EID"])

# At least 80% of encounter is ICU level-of-care -------------------------------
prop_ICU_time = rep(0, length(EIDs_sub))
for(i in 1:length(EIDs_sub)) {
    sub_loc = level_of_care_vec[data_format[,"EID"] == EIDs_sub[i]]
    prop_ICU_time[i] = mean(sub_loc == "Intensive Care")
}

EIDs_sub = EIDs_sub[prop_ICU_time > 0.8]
ind_icu_time = which(data_format[,"EID"] %in% EIDs_sub)
data_format = data_format[ind_icu_time, ]
level_of_care_vec = level_of_care_vec[ind_icu_time]
if(nrow(data_format) != length(level_of_care_vec)) print("error")

# Selection of training patients -----------------------------------------------
rbc_patients = unique(data_format[data_format[,"RBC_rule"] == 1,"EID"])

clinic_review = read.csv('../selected_encounters_update.csv')
clinical_id = as.numeric(clinic_review$key)
clinical_id_yes = clinical_id[clinic_review$Bleeding.event. %in% c("Yes", "Y")]
clinical_id_no = clinical_id[clinic_review$Bleeding.event. %in% c("No", "N")]

clinical_patients = clinical_id[clinical_id %in% EIDs_sub]
set.seed(2025)
TEST_SET = c(sample(clinical_patients[clinical_patients %in% clinical_id_no], 
                  size = 2, replace = F),
             sample(clinical_patients[clinical_patients %in% clinical_id_yes], 
                    size = 3, replace = F))
save(TEST_SET, file = "Data/TEST_SET.rda")

clinical_patients = clinical_patients[!(clinical_patients %in% TEST_SET)]
rbc_patients = rbc_patients[!(rbc_patients %in% TEST_SET)]
rbc_clinic_both = unique(c(rbc_patients, clinical_patients))

remaining_ids = EIDs_sub[!(EIDs_sub %in% c(rbc_clinic_both, TEST_SET))]

set.seed(2025)
EIDs_FINAL = c(rbc_clinic_both, sample(remaining_ids, size = 500 - length(rbc_clinic_both),
                      replace = F))
EIDs_FINAL = sort(EIDs_FINAL)

# double check TEST_SET is NOT in the training set!
if(sum(TEST_SET %in% EIDs_FINAL) != 0) print("ERROR: TEST IDs are in TRAIN IDs")

data_format = data_format[data_format[,"EID"] %in% EIDs_FINAL, ]
save(data_format, file = "Data/data_format_train.rda")

# ------------------------------------------------------------------------------
# FINAL FORMATTING OF MEDICATIONS ----------------------------------------------
# ------------------------------------------------------------------------------

# Finishing the medication formatting ------------------------------------------
load('Data/hr_cont_design.rda')
load('Data/hr_disc_design.rda')
load('Data/map_cont_design.rda')
load('Data/map_disc_design.rda')

med_ind_master = rep(0, length(EIDs_FINAL))
for(i in 1:length(EIDs_FINAL)) {
    med_ind_master[i] = which(EIDs == EIDs_FINAL[i])
}

hr_cont_design = hr_cont_design[med_ind_master]
hr_disc_design = hr_disc_design[med_ind_master]
map_cont_design = map_cont_design[med_ind_master]
map_disc_design = map_disc_design[med_ind_master]

# Dn_omega_names = c(hr_cont_names, hr_disc_names, map_cont_names, map_disc_names)
load('Data/Dn_omega_names.rda')
hr_cont_names = Dn_omega_names[1:ncol(hr_cont_design[[1]])]
start_ind = ncol(hr_cont_design[[1]])
hr_disc_names = Dn_omega_names[(start_ind + 1):(start_ind + ncol(hr_disc_design[[1]]))]
start_ind = ncol(hr_cont_design[[1]]) + ncol(hr_disc_design[[1]])
map_cont_names = Dn_omega_names[(start_ind + 1):(start_ind + ncol(map_cont_design[[1]]))]
start_ind = ncol(hr_cont_design[[1]]) + ncol(hr_disc_design[[1]]) + ncol(map_cont_design[[1]])
map_disc_names = Dn_omega_names[(start_ind + 1):(start_ind + ncol(map_disc_design[[1]]))]

# Matching the number of rows of the design matrices to the new data format ----
for(i in 1:length(EIDs_FINAL)) {
    sub_df = data_format[data_format[,"EID"] == EIDs_FINAL[i], ]
    
    if(!is.null(hr_cont_design[[i]])) {
        hr_cont_design[[i]] = hr_cont_design[[i]][1:nrow(sub_df), ,drop=F]
    } 
    
    if(!is.null(hr_disc_design[[i]])) {
        hr_disc_design[[i]] = hr_disc_design[[i]][1:nrow(sub_df), ,drop=F]
    } 
    
    if(!is.null(map_cont_design[[i]])) {
        map_cont_design[[i]] = map_cont_design[[i]][1:nrow(sub_df), ,drop=F]
    } 
    
    if(!is.null(map_disc_design[[i]])) {
        map_disc_design[[i]] = map_disc_design[[i]][1:nrow(sub_df), ,drop=F]
    } 
}

# Formatting Dn_omega ----------------------------------------------------------
Dn_omega = vector(mode = 'list', length = length(EIDs_FINAL))
hr_cont_num = length(hr_cont_names)
hr_disc_num = length(hr_disc_names)
map_cont_num = length(map_cont_names)
map_disc_num = length(map_disc_names)
total_cols_hr = hr_cont_num + hr_disc_num
total_cols_map = map_cont_num + map_disc_num
for(i in 1:length(Dn_omega)) {
    print(paste0(i, " in ", length(Dn_omega)))
    Dn_omega[[i]] = vector(mode = 'list', length = sum(data_format[,"EID"] == EIDs_FINAL[i]))
    for(j in 1:length(Dn_omega[[i]])) {
        med_mat = matrix(0, nrow = 4, ncol = total_cols_hr + total_cols_map)
        hr_info = map_info = NULL
        
        if(is.null(hr_cont_design[[i]])) {
            hr_info = rep(0, hr_cont_num)    
        } else {
            hr_info = hr_cont_design[[i]][j,]
        }
        
        if(is.null(hr_disc_design[[i]])) {
            hr_info = c(hr_info, rep(0, hr_disc_num))
        } else {
            hr_info = c(hr_info, hr_disc_design[[i]][j,])
        }
        names(hr_info) = NULL
        
        if(is.null(map_cont_design[[i]])) {
            map_info = rep(0, map_cont_num)    
        } else {
            map_info = map_cont_design[[i]][j,]
        }
        
        if(is.null(map_disc_design[[i]])) {
            map_info = c(map_info, rep(0, map_disc_num))
        } else {
            map_info = c(map_info, map_disc_design[[i]][j,])
        }
        names(map_info) = NULL
        
        med_mat[2, 1:total_cols_hr] = hr_info
        med_mat[3, (total_cols_hr+1):(total_cols_hr+total_cols_map)] = map_info
        
        Dn_omega[[i]][[j]] = med_mat
    }
}

save(Dn_omega, file = paste0('Data/Dn_omega.rda'))

# Understanding what the mean of Dn_omega should be
load("Data/med_select_FINAL.rda")
upp_down_omega = matrix(nrow = length(Dn_omega_names), ncol = 2)
upp_down_omega[,1] = Dn_omega_names
ind = 1
for(i in 1:length(hr_cont_names)) {
    if(upp_down_omega[ind,1] != hr_cont_names[i]) {
        print(hr_cont_names[i])
    } else{
        ef = unique(med_select_FINAL$hr[med_select_FINAL$med_name_admin == hr_cont_names[i]]) 
        print(ef)
        upp_down_omega[ind, 2] = as.numeric(ef)
    }
    ind = ind + 1
}
for(i in 1:length(hr_disc_names)) {
    if(upp_down_omega[ind,1] != hr_disc_names[i]) {
        print(hr_disc_names[i])
    } else{
        ef = unique(med_select_FINAL$hr[med_select_FINAL$med_name_admin == hr_disc_names[i]]) 
        print(ef)
        upp_down_omega[ind, 2] = as.numeric(ef)
    }
    ind = ind + 1
}
for(i in 1:length(map_cont_names)) {
    if(upp_down_omega[ind,1] != map_cont_names[i]) {
        print(map_cont_names[i])
    } else{
        ef = unique(med_select_FINAL$map[med_select_FINAL$med_name_admin == map_cont_names[i]]) 
        print(ef)
        upp_down_omega[ind, 2] = as.numeric(ef)
    }
    ind = ind + 1
}
for(i in 1:length(map_disc_names)) {
    if(upp_down_omega[ind,1] != map_disc_names[i]) {
        print(map_disc_names[i])
    } else{
        ef = unique(med_select_FINAL$map[med_select_FINAL$med_name_admin == map_disc_names[i]]) 
        print(ef)
        upp_down_omega[ind, 2] = as.numeric(ef)
    }
    ind = ind + 1
}
rm(med_select_FINAL)
mean_dn_omega = as.numeric(c(upp_down_omega[,2]))
print(c(mean_dn_omega))
cat(mean_dn_omega, sep = ',')
# -1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
# -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
# -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
# -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1