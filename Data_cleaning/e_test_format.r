load('Data/data_format.rda')
load('Data/level_of_care.rda')
load('Data/cov_info.rda')  
load('Data/TEST_SET_update.rda')
TEST_SET = sort(TEST_SET)

EIDs = unique(data_format[,"EID"])

# ------------------------------------------------------------------------------
# FINAL FORMATTING OF TEST SET -------------------------------------------------
# ------------------------------------------------------------------------------

ind_dat = which(data_format[,"EID"] %in% TEST_SET)
data_format = data_format[ind_dat, ]
level_of_care_vec = level_of_care_vec[ind_dat]
if(nrow(data_format) != length(level_of_care_vec)) print("error")

# Adding RBC rule (if >= 3 RBCs in 12hr window) --------------------------------

min_three_RBC = unique(data_format[data_format[,"n_RBC_admin"] >= 3, "EID"])
bleed_pat = min_three_RBC

for(i in 1:length(bleed_pat)) {
    sub_dat = data_format[data_format[,"EID"] == bleed_pat[i], ]
    
    # Check in any 12 hour period
    max_time = tail(sub_dat[,"time"], 1)
    when_rbc = c(1, which(diff(sub_dat[,"n_RBC_admin"]) != 0))
    
    for(j in 1:length(when_rbc)) {
        s_time = sub_dat[when_rbc[j], "time"]
        e_time_12 = s_time + 720
        RBC_diff_12 = 0
        
        if (e_time_12 <= max_time) {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            ind_12 = order(abs(sub_dat[,"time"] - e_time_12))[1]
            RBC_diff_12 = sub_dat[ind_12, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        } else {
            s_ind = order(abs(sub_dat[,"time"] - s_time))[1]
            e_ind = order(abs(sub_dat[,"time"] - max_time))[1]
            RBC_diff_12 = sub_dat[e_ind, "n_RBC_admin"] - sub_dat[s_ind, "n_RBC_admin"]
        }
        
        if(RBC_diff_12 >=3) {
            data_format[data_format[,"EID"] == bleed_pat[i], "RBC_rule"] = 1
            break
        }
    }
}

save(data_format, file = 'Data/data_format_TEST.rda')
# Finishing the medication formatting ------------------------------------------

load('Data/hr_cont_design.rda')
load('Data/hr_disc_design.rda')
load('Data/map_cont_design.rda')
load('Data/map_disc_design.rda')

med_ind_master = rep(0, length(TEST_SET))
for(i in 1:length(TEST_SET)) {
    med_ind_master[i] = which(EIDs == TEST_SET[i])
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
for(i in 1:length(TEST_SET)) {
    sub_df = data_format[data_format[,"EID"] == TEST_SET[i], ]
    
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
Dn_omega = vector(mode = 'list', length = length(TEST_SET))
hr_cont_num = length(hr_cont_names)
hr_disc_num = length(hr_disc_names)
map_cont_num = length(map_cont_names)
map_disc_num = length(map_disc_names)
total_cols_hr = hr_cont_num + hr_disc_num
total_cols_map = map_cont_num + map_disc_num
for(i in 1:length(Dn_omega)) {
    print(paste0(i, " in ", length(Dn_omega)))
    Dn_omega[[i]] = vector(mode = 'list', length = sum(data_format[,"EID"] == TEST_SET[i]))
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

save(Dn_omega, file = paste0('Data/Dn_omega_TEST.rda'))

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

# 500 subject "Final count of all training data = 1516"
# "Number of possible subjects: 1441"
# -1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
# -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
# -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
# -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1

# 1000 subject "Final count of all training data = 1516"
# "Number of possible subjects: 1441"
# -1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
# -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
# -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
# -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1

# Test Set
# -1, 1, 1,-1,-1, 1, 1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,-1,-1,-1, 1,
# -1, 1,-1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1, 1, 1, 1,-1, 1,
# -1,-1,-1, 1,-1, 1,-1, 1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,
# -1,-1, 1, 1, 1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1

