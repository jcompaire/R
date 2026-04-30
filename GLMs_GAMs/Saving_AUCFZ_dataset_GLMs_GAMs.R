# Clear console and environment ####
cat("\014") 
rm(list=ls()) 
graphics.off()
gc()
# Load packages ####
pk <- c("comprehenr", "lubridate", "ncdf4", "readxl", "tidyverse")
lapply(pk, require, character.only = TRUE)
# Set working directory ####
cwd <- getwd()
if(grepl('Users/Worm',cwd, fixed = TRUE)){
  setwd("~/GoogleDrive/CIMA/CTMFM/")
  nc_wd = paste0(getwd(),"datasets_merged_ENVIRONMENT_FISH/")
    # (
    # "/Users/Worm/GoogleDrive/CIMA/CTMFM/datasets_merged_ENVIRONMENT_FISH/"
  # )
  # path = '/Users/Worm/GoogleDrive/CIMA/CTMFM/2_GLMs_GAMs/GLMs/'
} else{
  setwd("~/gdrive/CIMA/CTMFM/")
  # nc_wd = (
  #   "/home/gsus/gdrive/CIMA/CTMFM/datasets_merged_ENVIRONMENT_FISH/"
  # )
  # path = '/home/gsus/gdrive/CIMA/CTMFM/2_GLMs_GAMs/GLMs/'
}
nc_wd = paste0(getwd(),"/datasets_merged_ENVIRONMENT_FISH/")
path = paste0(getwd(),"/2_GLMs_GAMs/")
# === AUCFZ POLYGON COORDINATES ====
# Specify sheet by number
aucfz <- read_excel("files_EXCEL/AUCFZ_coordinates.xlsx", sheet = 1)
# === SPAWNING POLYGON COORDINATES ====
# Specify sheet by number
spawning <- read_excel("files_EXCEL/Mfurnieri_spawning_polygon.xlsx", sheet = 1)
# === LANDINGS and ENVIRONMENTAL data LAGGED ====
data_nc <- nc_open(paste0(
  nc_wd, 'Mfurnieri_lagged_10yr_DATAFRAME_raw_continuous_standardized.nc'))
names_sp <- (names(data_nc$var)) # Species' names
names_col <- (data_nc[["dim"]][["variables"]][["vals"]]) # Env + Fish names
# names_col <- names_col[-grep("chl",names_col)]
no_chl <- names_col[-grep("chl",names_col)]
idx_no_chl <- to_vec(for(i in 1:length(names_col))
  str_which(names_col, no_chl[i], negate = FALSE))
names_lag <- (data_nc[["dim"]][["lags"]][["vals"]]) # Lag's names
sp <- ncvar_get(data_nc, names_sp[1]) # Mfurnieri
sp <- sp[idx_no_chl, , ]
names_col <- names_col[idx_no_chl] # names_col without chl variables
# 
names_time <- as.Date(
  (data_nc[["dim"]][["time"]][["vals"]]), origin = '1986-01-01')
# 
data.set <- list()
data.set[[1]] <- aucfz
data.set[[2]] <- spawning
data.set[[3]] <- sp
data.set[[4]] <- names_col
data.set[[5]] <- names_time
data.set[[6]] <- names_lag

# Save dateset
save(data.set, file="glm_gam_aucfz_data.RData")
