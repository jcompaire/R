# Clear console and environment ####
cat("\014") 
rm(list=ls()) 
# Load packages needed for analysis and figures ####
pk <- c("comprehenr", "lubridate", "ncdf4", "tidyverse", "WaveletComp")
lapply(pk, require, character.only = TRUE)
# Set working directory, invoke functions and load datasets  ####
setwd("/Users/Worm/GoogleDrive/CIMA/Environmental_data/ncfiles_fish_env/") 
# 
nc_data <- nc_open("Fishes_lagged_120months_DATAFRAME_raw.nc")
time <- ncvar_get(nc_data, "time")
tunits <- ncatt_get(nc_data, "time", "units") #check units
dim(time) # check dimension
time_obs <- as.POSIXct(time*86400, origin = "1986-01-01 00:00:00",
                       format = c("%Y-%m-%d %H:%M:%OS"), tz = 'GMT')
dim(time_obs) #should be 408
range(time_obs)
# fillvalue <- ncatt_get(nc_data, "C.guatucupa", "_FillValue")
# C.guatucupa ####
vrs <- nc_data$var$Cynoscion_guatucupa[["dim"]][[1]][["vals"]]
print(vrs)
f_array <- ncvar_get(nc_data, "Cynoscion_guatucupa")
f_array <- aperm(f_array, c(2,3,1)) # reordering array
landings <- f_array[, 1, 22] # landings
pf <- plot(ts(landings, start = c(1986,1), frequency = 12))
# Lags with more significant correlations: lag 2 RIVER, lag 10 CHL
# whatever lag + 1 
riv <- f_array[ ,3 , 5]
chl <- f_array[ ,11 , 17]
df <- as.data.frame(cbind(landings, riv, chl))
df[["date"]] <- as.POSIXct(time_obs)
df$Common_name <- "Stripped weakfish"
df$Scientific_name <- "C. guatucupa"
df <- df[c(4, 6, 5, 1:3)]
cguatucupa <- df[complete.cases(df[,c("landings")]),]
# Clear console and environment keeping dataset, packages and functions ####
cat("\014") 
rm(list = setdiff(ls(), (c("nc_data", "time_obs", "pk", lsf.str(),
                           "cguatucupa"))))
gc() # Free unused memory
# M.furnieri ####
vrs <- nc_data$var$Micropogonias_furnieri[["dim"]][[1]][["vals"]]
print(vrs)
f_array <- ncvar_get(nc_data, "Micropogonias_furnieri")
f_array <- aperm(f_array, c(2,3,1)) # reordering array
landings <- f_array[, 1, 22] # landings
pf <- plot(ts(landings, start = c(1986,1), frequency = 12))
# Lags with more significant correlations: lag 8 KD, lag 3 SSS, lag 14 WMI
# whatever lag + 1 
kd  <- f_array[ ,9 , 11]
sss <- f_array[ ,4 , 14]
wmi <- f_array[ ,15 , 20]
df <- as.data.frame(cbind(landings, kd, sss, wmi))
df[["date"]] <- as.POSIXct(time_obs)
df$Common_name <- "Whitemouth croaker"
df$Scientific_name <- "M. furnieri"
df <- df[c(5, 7, 6, 1:4)]
mfurnieri <- df[complete.cases(df[,c("landings")]),]
# Clear console and environment keeping dataset, packages and functions ####
cat("\014") 
rm(list = setdiff(ls(), (c("nc_data", "time_obs", "pk", lsf.str(),
                           "cguatucupa", "mfurnieri"))))
gc() # Free unused memory
# M.hubbsi ####
vrs <- nc_data$var$Merluccius_hubbsi[["dim"]][[1]][["vals"]]
print(vrs)
f_array <- ncvar_get(nc_data, "Merluccius_hubbsi")
f_array <- aperm(f_array, c(2,3,1)) # reordering array
landings <- f_array[, 1, 22]
pf <- plot(ts(landings, start = c(1986,1), frequency = 12))
# Lags with more significant correlations: lag 5 RIVER, lag 0 CHL, lag 11 WMI
# whatever lag + 1  
riv <- f_array[ ,6 , 5]
chl <- f_array[ ,1 , 17]
wmi <- f_array[ ,12 , 20]
df <- as.data.frame(cbind(landings, riv, chl, wmi))
df[["date"]] <- as.POSIXct(time_obs)
df$Common_name <- "Argentine hake"
df$Scientific_name <- "M. hubbsi"
df <- df[c(5, 7, 6, 1:4)]
mhubbsi <- df[complete.cases(df[,c("landings")]),]
# Average wavelet power for environmental variables ####
vrs <- nc_data$var$Cynoscion_guatucupa[["dim"]][[1]][["vals"]]
print(vrs)
f_array <- ncvar_get(nc_data, "Cynoscion_guatucupa")
f_array <- aperm(f_array, c(2,3,1)) # reordering array
lg <- 1 # first element 
riv <- f_array[ ,lg , 5]
kd  <- f_array[ ,lg , 11] 
sss <- f_array[ ,lg , 14] 
chl <- f_array[ ,lg , 17] 
wmi <- f_array[ ,lg , 20] 
df <- as.data.frame(cbind(riv, kd, sss, chl, wmi))
df[["date"]] <- as.POSIXct(time_obs)
environmental <- df[c(6, 1:5)]
# Save query species one by one and environmental dataset in one object ####
data.set <- list()
data.set[[1]] <- cguatucupa
data.set[[2]] <- mfurnieri
data.set[[3]] <- mhubbsi
data.set[[4]] <- environmental
# Save raw dateset
save(data.set, file="wavelet_data.RData")
