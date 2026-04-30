library(diffr)
getwd()
setwd("/Users/Worm/GoogleDrive/CIMA/CTMFM/GLMs/")

fil1 = "/Users/Worm/Downloads/time_series_aucfz_functions.R"
fil2 = "/Users/Worm/GoogleDrive/GitHub/R/TimeSeries/time_series_aucfz_functions.R"

diffr(fil1, fil2, before = "f1", after = "f2")

fil1 <- '~/gdrive/CIMA/CTMFM/3_Clustering_TimeSeries/GLMs_proofs.R'
fil2 <- '~/gdrive/CIMA/CTMFM/scripts_py/try_glm_proof.R'
diffr(fil1, fil2, before = "f1", after = "f2")
