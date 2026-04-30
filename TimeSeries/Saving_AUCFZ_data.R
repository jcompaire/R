# Clear console and environment ####
cat("\014"); rm(list=ls())
# Load packages ####
pk <- c("comprehenr", "lubridate", "ncdf4", "readxl", "tidyverse")
lapply(pk, require, character.only = TRUE)
# Set working directory ####
setwd("~/gdrive/CIMA/CTMFM/")
setwd("~/GoogleDrive/CIMA/CTMFM/")
# === AUCFZ POLYGON COORDINATES ====
# Specify sheet by number
aucfz <- read_excel("AUCFZ_coordinates.xlsx", sheet = 1)
# === CPUE vs COMMERCIAL LANDINGS ====
# Specify sheet by number
df_cl <- read_excel("DINARA_Landings_CPUE.xlsx", sheet = 2)
# === LANDINGS =====
# Specify sheet by name ####
# data.set <- read_excel("", sheet = "NAME")
# Specify sheet by number
df1 <- read_excel("CTMFM_Landings_1996_2019_sorted_dates.xlsx", sheet = 1)
df2 <- read_excel("CTMFM_Landings_2020_2021_sorted_dates.xlsx", sheet = 1)
df <- rbind(df1,df2) # Merge both dataset by row
rm(df1,df2)
# To use dates, these should be provided in numeric format: yyyy/mm/dd
month_n <- match(df$Month, month.abb) # get month number column
df$date_num <- ymd(paste(df$Year,month_n, 1 , sep="-")) # Year-Month-Day
rm(month_n)
df <- df[,c(8,2,3,4,5,6,7)] # reorder columns
# df <- df[,c(1,2,3,8,4,5,6,7)] # reorder columns
# Changing common names from Spanish to English
df$Common_name[df$Common_name == 'Corvina'] <- 'Whitemouth croaker'
df$Common_name[df$Common_name == 'Pescadilla'] <- 'Stripped weakfish'
df$Common_name[df$Common_name == 'Merluza'] <- 'Argentine hake'
# Crop scientific names
df$Scientific_name[df$Scientific_name == 'Micropogonias furnieri'] <- 'M. furnieri'
df$Scientific_name[df$Scientific_name == 'Cynoscion guatucupa'] <- 'C. guatucupa'
df$Scientific_name[df$Scientific_name == 'Merluccius hubbsi'] <- 'M. hubbsi'
# Get species names to use dataset for each SPECIES ###
sp_names <- df %>% select(5)
ux <- unique(sp_names); rm(sp_names)
sname <- c('C. guatucupa','M. furnieri','M. hubbsi')
## Subset by country ####
Uru <- subset(df, Country == "Uruguay")
Arg <- subset(df, Country == "Argentina")
# Save query species one by one in a whole ordered dataset ####
data.set = list(aucfz, df_cl)
for (j in 1:length(sname)) {
  # Select query specie
  Q_sp <- which(toString(sname[j]) == Uru[ ,5])
  Q_spA <- which(toString(sname[j]) == Arg[ ,5])
  # Save to GENERAL SUBPLOT for EACH SPECIES
  data.set[[j+2]] <- Uru[Q_sp,]
  data.set[[j+5]] <- Arg[Q_spA,]
}
# Save raw dateset
# save(data.set, file="time_series_aucfz_data.RData")
# Check series
# t <- data.set[[4]]
# tt <-(ts(t[c("Landing")],
#               start = c(1996,1), end = c(2021,12), frequency = 12))
# plot(decompose(tt))
# Clear console and environment keeping dataset, packages and functions ####
cat("\014") 
rm(list = setdiff(ls(), (c("data.set", "pk", lsf.str()))))
gc() # Free unused memory
# === HEATMAP DATA ====
# Specify sheet by number
df <- read_excel("HEATMAP_values.xlsx", sheet = 1)
df$species <- as.factor(df$species)
df$variable <- as.factor(df$variable)
j <- length(data.set)
data.set[[j+1]] <- df
# Clear console and environment keeping dataset, packages and functions ####
cat("\014") 
rm(list = setdiff(ls(), (c("data.set", "pk", lsf.str()))))
gc() # Free unused memory
# === WAVELET DATA ====
# # Load packages needed for analysis and figures ####
# pk <- c("comprehenr", "lubridate", "ncdf4", "tidyverse", "WaveletComp")
# lapply(pk, require, character.only = TRUE)
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
rm(list = setdiff(ls(), (c("data.set", "nc_data", "time_obs", "pk", lsf.str(),
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
rm(list = setdiff(ls(), (c("data.set","nc_data", "time_obs", "pk", lsf.str(),
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
j <- length(data.set)
data.set[[j+1]] <- cguatucupa
data.set[[j+2]] <- mfurnieri
data.set[[j+3]] <- mhubbsi
data.set[[j+4]] <- environmental
# Save dateset
# save(data.set, file="wavelet_data.RData")
save(data.set, file="time_series_wavelet_aucfz_data2.RData")
