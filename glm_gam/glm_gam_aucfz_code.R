## Code description ####
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Author: Jesus C. Compaire
## Institution: Centro de Investigaciones del Mar y la Atmósfera (CIMA)
## Position: Postdoctoral researcher
## Contact details: jesus.canocompaire@uca.es
## Date created: Dec-2023
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Code to replicate the statistical analyses and generate figures 
## performed in the manuscript:
## Compaire, J.C., Simionato C.G.  Moreira, D. & Acha, E.M. (2023).
## Modeling environmental effects on fishery landings: A case study of 
## whitemouth croaker (Micropogonias furnieri) in the Southwestern
## Atlantic shelf
##
## -- -- -- -- -- -- -- -- -- -- -- -- --
#
# Clear console, environment, figures and free memory ####
cat("\014"); rm(list=ls()); graphics.off(); gc()
# Load packages needed for analysis and figures ####
pk <- c("comprehenr", "dplyr", "effects", "formula.tools", "ggplot2", "ggpubr",
        "ggthemes", "ggspatial", "gratia",  "grid", "gridExtra", "hrbrthemes",
        "mgcv", "MuMIn", "ncdf4", "patchwork", "performance", "plyr",
        "readxl","rnaturalearth", "rlist", "scales", "stats", "stringr")
lapply(pk, require, character.only = TRUE)
# Set working directory, invoke functions and load datasets ####
setwd("~/gdrive/GitHub/R/GLMs_GAMs/") # directory
source("glm_gam_aucfz_functions.R") # functions
load("glm_gam_aucfz_data.RData") # data
output <- paste0(getwd(), "/output_figures/") # Directory to put figures
# ============ MAP ======================================================= ####
aucfz <- as.data.frame(data.set[1]) # AUCFZ polygon
spawning <- as.data.frame(data.set[2]) # Spawning ground polygon
rs <- raster("elevation_AUCFZ_ETOPO2022.tiff") # Bathymetry
## Preparing raster object to plot with ggplot2
my_region <- extent(-60, -50, -40, -33.1)
r_points <- rasterToPoints(crop(rs, my_region)) # from raster to points
#  The error that appears at the first time is an "error message" but not a real
#  error in the strict sense as the script continues to run and the results are
#  are not affected. For details, see: https://stackoverflow.com/a/61623381
#
bathy <- data.frame(r_points); rm(r_points) # convert to dataframe
colnames(bathy) <- c("lon", "lat", "elevation")
# Breaks and labels for color scale
zmin <- min(bathy$elevation)
brks <- c(0, -50, -200, -1000, -3000, -5000, zmin)
clabs <- c("(0-50]", "(50-200]", "(200-1000]", "(1000-3000]",
           "(3000-5000]", "> 5000")
bathy$discrete <- cut(bathy$elevation, breaks = brks, labels = clabs) # discrete
head(bathy)
# Saving plot
file_name = paste0(file = output , "Fig1_Map",
                   ".pdf", sep="")
pdf(file_name, width=9, height=6, bg = "transparent", colormodel = "cmyk",
    compress = T)
par(mfcol=c(1,1))
p1 <- map_aucfz(df = bathy, plygn = aucfz, plygn2 = spawning,
                scale = "discrete", clabs)
# South America map
sa_countries <- ne_countries(continent = "South America", returnclass = "sf")
p2 <- ggplot() +
  geom_sf(data=sa_countries, color = "grey22", fill = "grey90") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size=0.5),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(-0.3, -0.3, -0.3, -0.3, unit = "cm")) +
  # Bounding box
  geom_rect(aes(xmin = -60, xmax = -50, ymin = -40, ymax = -33),
            color = "black", lwd = 0.6, fill = NA) +
  ## Text
  annotate(
    geom = "text", x = -53, y = -7, label = "South America", 
    fontface = "plain", family = "Times", color = "grey22", size = 2.5) +
  annotate(
    geom = "text", x = -50, y = -48, label = "Atlantic Ocean",
    fontface = "plain", family = "Times", color = "grey22", size = 2.5)
# Create inset map
grid.newpage()
v1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # larger map
v2 <- viewport(width = 0.23, height = 0.23, x = 0.915, y = 0.85)  # inset map
print(p1, vp = v1)
print(p2, vp = v2)
dev.off()
# Clear console and environment keeping data, packages and functions ####
cat("\014")
ds <- c("path", "data.set", "output")
rm(list = setdiff(ls(), (c(ds, "pk", lsf.str()))))
gc() # Free unused memory
## Getting values for GLM and GAM analyses =================== ####
sp <- data.set[[3]]
names_col <- data.set[[4]]
names_dates <- data.set[[5]]
names_lag <- data.set[[6]]
# ============ GLMs ANALYSES ============================================= ####
## GLMs analyses: TH covariates - TOTAL landings ####
# Diverting the output into a single external file
sink(paste0(file = 'GLMs_Mfurnieri_outputs.txt')) 
# Evaluating different families (fmly) and percentiles (th)
fmly <- list()
fmly[[1]] <- Gamma(link = "inverse")
fmly[[2]] <- inverse.gaussian(link = "1/mu^2")
th <- c('25','50','75')
# i) === Results for each family ===
glms_summary_fm <- list()
for(f in 1:length(fmly)){
  # print(paste0(fmly[[f]][[1]],'(link =',fmly[[f]][[2]], ')'))
  # ii) === Results for each percentile ===
  glms_summary_th <- list()  
  for(k in 1:length(th)){
    # print(th[k])
    # List comprehension to get index (columns number) of variables
    vq <- c(th[k], 'total')
    idx_var <- to_vec(for(i in 1:length(vq))
      str_which(names_col, vq[i], negate = FALSE))
    names_var <- names_col[idx_var] # variable names complete
    # Getting climatic indices by pairs
    ci <- c('soi','sam','enso')
    cixs <- to_list(for(i in 1:length(ci)) 
      c(combn(ci, i, simplify = FALSE)))
    cixss <- list()
    cixss[[1]] <- cixs[[2]] # pairs of climatic indexes
    cixss[[1]][[4]] <- unlist(cixs[[3]]) # all climatic indexes
    names_var_ix <- NA # variable names short for each climatic index
    # iii) === Results for each climatic index and without them ===
    glms_summary_ix <- list() 
    for(ix in 1:length(cixss[[1]])){
      # Getting two of the climatic indices and 'th' landings to omit them
      cith <- paste(paste(cixss[[1]][[ix]], collapse = '|'), 'th', sep = '|')
      names_var_ix = names_var[
        -grep(pattern = cith, names_var)]
      # List comprehension to get index (columns number) of variables
      idx_ss <- to_vec(for(i in 1:length(names_var_ix))
        str_which(names_col, names_var_ix[i], negate = FALSE))
      # Removing tag "th"
      tag <- c(th[k])
      for(i in 1:length(tag)){
        names_var_ix <- gsub(tag[i],'', names_var_ix)
      }
      # print(names_var_ix[1])
      # Covariates and response variable names
      xnames <- head(names_var_ix, -1); ynames <- tail(names_var_ix, 1)
      # GLMs formulas ADDITIVE glms_f <- comb_glm(xnames, ynames) 
      # GLMs formulas INTERACTION
      inter <- c("river", "kd", "sss") # interaction terms
      glms_f <- comb_glmi(xnames, ynames, inter) # getting equations
      if(
        ix < length(cixss[[1]])
      ){
        glms_f <- glms_f[grep(
          pattern = paste(cixss[[1]][[4]], collapse = '|'), (glms_f))]
      } else {
        glms_f <- glms_f
      }
      #  To divert the output into separate external files ordered by family,
      #  percentile and climatic index or river:
      #  - Comment lines 74, 172 and 173
      #  - Uncomment lines 135-138 and 164-165
      # fmly_name <- paste0(fmly[[f]][[1]])
      # sink()
      # sink(paste0(file = path, 'GLMs_Mfurnieri_', fmly_name, '_',
      #             th[k], 'th_', toupper(names_var_ix[1]), '.txt'))
      # iv) === Results for each lag ===
      glms_summary_lg <- list() 
      for(i in 1:length(names_lag)){
        spq <- t(sp[c(idx_ss), , i]) # dataset to each lag
        spq <- as.data.frame(spq)
        # spq$V1 <- as.factor(spq$V1) # climatic index as qualitative variable
        colnames(spq) <- names_var_ix # adding column names to DF
        # v) === Results for each equation ===
        glms_summary_eq <- matrix(NA, ncol = 4, nrow = length(glms_f))
        for(j in 1:length(glms_f)){
          eq <- glms_f[[j]]
          output <- try_glm(eq, spq, i-1, j, fmly[[f]], th[k])
          glms_summary_eq[j,] <- c(output$lag, output$eq ,output$DE, output$AICc)
          colnames(glms_summary_eq) <- c('lag','eq','DE','AICc')
          glms_summary_eq <- data.frame(glms_summary_eq)
          glms_summary_eq <- na.omit(glms_summary_eq)
        }
        # Calculate ∆AICc only whether there are valid models
        if (dim(glms_summary_eq)[1] != 0)
        {
          glms_summary_eq$AICc <- (glms_summary_eq$AICc -
                                     min(glms_summary_eq$AICc))
        }
        glms_summary_lg <- rbind(glms_summary_lg,glms_summary_eq)
      }
      # sink() # sets the output back to the console
      # closeAllConnections() # to print on the console again
      glms_summary_ix[[ix]] <- glms_summary_lg
    }
    glms_summary_ix <- setNames(glms_summary_ix, c(rev(ci), 'river'))
    glms_summary_th[[k]] <- glms_summary_ix
  }
  glms_summary_th <- setNames(glms_summary_th, c(th))
  glms_summary_fm[[f]] <- glms_summary_th 
}
sink() # sets the output back to the console
closeAllConnections() # to print on the console again
glms_summary_fm <- setNames(glms_summary_fm, 
                            c(fmly[[1]]$family,fmly[[2]]$family))
save(glms_summary_fm, file = paste0('GLMs_summary_fm.RData'))
## GLMs plot summarized ####
load(paste0("GLMs_summary_fm.RData"))
# Checking percentiles from output file (.txt) to discard those models with
# high correlation among variables (VIF > 10)
ix <- c('enso','sam','soi','river')
# i)  50th --------------------------------------------------------------------
th <- "50"
# === Gamma ===
# Visualize number of GLMs according to each index and their DE values
to_vec(for(i in 1:length(ix))
  sapply(glms_summary_fm[["Gamma"]][[th]][[ix[i]]][1], NROW))
to_vec(for(i in 1:length(ix))
  (glms_summary_fm[["Gamma"]][[th]][[ix[i]]][['DE']]))
# All VIF values  < 10
lag_v <- vector_de_lag(glms_summary_fm, 'Gamma',
                       th, ix, 'lag')
de_v <- vector_de_lag(glms_summary_fm, 'Gamma',
                      th, ix, 'DE')
t50g <- table_de_lag(de_v, lag_v, 'gamma') 
# === Inverse gaussian ===
to_vec(for(i in 1:length(ix))
  sapply(glms_summary_fm[["inverse.gaussian"]][[th]][[ix[i]]][1], NROW))
to_vec(for(i in 1:length(ix))
  (glms_summary_fm[["inverse.gaussian"]][[th]][[ix[i]]][['DE']]))
# All VIF values  < 10
lag_v <- vector_de_lag(glms_summary_fm, 'inverse.gaussian',
                      th, ix, 'lag')
de_v <- vector_de_lag(glms_summary_fm, 'inverse.gaussian',
                      th, ix, 'DE')
t50i <- table_de_lag(de_v, lag_v, 'inverse.gaussian')
# Merging them into a single DF
cols_num <- c('max','avg','sd')
t50g[cols_num] <- sapply(t50g[cols_num], as.numeric)
t50i[cols_num] <- sapply(t50i[cols_num], as.numeric)
th50 <- rbind (t50g, t50i)
th50$avg <- ifelse(th50$avg == th50$max, NA, th50$avg)
# ii)  25th -------------------------------------------------------------------
th <- "25"
# === Gamma ===
# Visualize number of GLMs according to each index and their DE values
to_vec(for(i in 1:length(ix))
  sapply(glms_summary_fm[["Gamma"]][[th]][[ix[i]]][1], NROW))
to_vec(for(i in 1:length(ix))
  (glms_summary_fm[["Gamma"]][[th]][[ix[i]]][['DE']]))
# VIF > 10: rows 1, 3, 4 have to be dropped
rm_rws <- c(1, 3, 4)
lag_v <- vector_de_lag(glms_summary_fm, 'Gamma',
                       th, ix, 'lag')[-c(rm_rws)]
de_v <- vector_de_lag(glms_summary_fm, 'Gamma',
                      th, ix, 'DE')[-c(rm_rws)]
t25g <- table_de_lag(de_v, lag_v, 'gamma')  
# === Inverse gaussian ===
# Visualize number of GLMs according to each index and their DE values
to_vec(for(i in 1:length(ix))
  sapply(glms_summary_fm[['inverse.gaussian']][[th]][[ix[i]]][1], NROW))
to_vec(for(i in 1:length(ix))
  (glms_summary_fm[["inverse.gaussian"]][[th]][[ix[i]]][['DE']]))
vector_de_lag(glms_summary_fm, 'inverse.gaussian', th,  ix, 'DE')
# VIF > 10: rows 1, 2, 3 have to be dropped
rm_rws <- c(1, 2, 3)
lag_v <-vector_de_lag(glms_summary_fm, 'inverse.gaussian', 
                      th, ix, 'lag')[-c(rm_rws)]
de_v <- vector_de_lag(glms_summary_fm, 'inverse.gaussian',
                      th, ix, 'DE')[-c(rm_rws)]
t25i <- table_de_lag(de_v, lag_v, 'inverse.gaussian')
# Merging them into a single DF
t25g[cols_num] <- sapply(t25g[cols_num], as.numeric)
t25i[cols_num] <- sapply(t25i[cols_num], as.numeric)
th25 <- rbind(t25g, t25i)
th25$avg <- ifelse(th25$avg == th25$max, NA, th25$avg)
# iii)  75th ------------------------------------------------------------------
th <- "75"
# === Gamma ===
# Visualize number of GLMs according to each index and their DE values
to_vec(for(i in 1:length(ix))
  sapply(glms_summary_fm[["Gamma"]][[th]][[ix[i]]][1], NROW))
to_vec(for(i in 1:length(ix))
  (glms_summary_fm[["Gamma"]][[th]][[ix[i]]][['DE']]))
# VIF values  < 10
lag_v <- vector_de_lag(glms_summary_fm, 'Gamma',
                       th, ix, 'lag')
de_v <- vector_de_lag(glms_summary_fm, 'Gamma',
                      th, ix, 'DE')
t75g <- table_de_lag(de_v, lag_v, 'gamma')
# === Inverse gaussian ===
to_vec(for(i in 1:length(ix))
  sapply(glms_summary_fm[['inverse.gaussian']][[th]][[ix[i]]][1], NROW))
to_vec(for(i in 1:length(ix))
  (glms_summary_fm[["inverse.gaussian"]][[th]][[ix[i]]][['DE']]))
# All VIF values  < 10 VIF < 10
lag_v <-vector_de_lag(glms_summary_fm, 'inverse.gaussian',
                      th, ix, 'lag')
de_v <- vector_de_lag(glms_summary_fm, 'inverse.gaussian',
                      th, ix, 'DE')
t75i <- table_de_lag(de_v, lag_v, 'inverse.gaussian')
# Merging them into a single DF
t75g[cols_num] <- sapply(t75g[cols_num], as.numeric)
t75i[cols_num] <- sapply(t75i[cols_num], as.numeric)
th75 <- rbind (t75g, t75i)
th75$avg <- ifelse(th75$avg == th75$max, NA, th75$avg)
# Plot ------------------------------------------------------------------------
p50 <- plot_de_lag(th50) +
  theme(
    legend.direction = "vertical",
    legend.position = c(0.25,0.8),
    legend.key.size = unit(0.75, "cm"),
    legend.text = element_text(size =  12)
    # legend.title = element_text(size = 15, face = "bold")
  ) +
  annotate(geom = "text", x = 2, y = 90, label = "50th",
           color = "black", size = 10, family = "Times") +
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.12, 0.975),
    plot.tag = element_text(
      face = "bold", size = 20)
  ) 
p25 <- plot_de_lag(th25) +
  annotate(geom = "text", x = 2, y = 90, label = "25th",
           color = "black", size = 10, family = "Times") +
  labs(tag = "B") +
  theme(
    plot.tag.position = c(0.12, 0.975),
    plot.tag = element_text(
      face = "bold", size = 20)
  )
p75 <- plot_de_lag(th75) +
  annotate(geom = "text", x = 2, y = 90, label = "75th",
           color = "black", size = 10, family = "Times") +
  labs(tag = "C") +
  theme(
    plot.tag.position = c(0.12, 0.975),
    plot.tag = element_text(
      face = "bold", size = 20)
  )
# Merging plots ---------------------------------------------------------------
file_name = paste0(file = output , "Fig2_GLMs_summarized",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
figure <- ggarrange(p50, p25, p75, ncol = 3, nrow = 1, common.legend = F)
annotate_figure(
  figure, 
  bottom = textGrob(
    expression(paste("Lag (years)")),
    rot = 0, vjust = 0.3,
    gp = gpar(cex = 1.5, fontfamily = "Times" , fontsize = 12)
  ),
  left = textGrob(
    expression(paste("Deviance explained (%)")),
    rot = 90, vjust = 1.3,
    gp = gpar(cex = 1.5, fontfamily = "Times" , fontsize = 12)
  )) 
dev.off()
#
## Checking some models manually ####
# k = 25th: 1 - 50th: 2 - 75th: 3  / ix = ENSO: 1 - SAM: 2 - SOI: 3 - RIVER: 4
# l = lag (from 0 to 10) / j = number of equation (from "glms_summary_fm")
th <- c('25','50','75')
# i.1) Lag 1 - 50th - (river + kd + river:kd) ---------------------------------
hf <- handyf_glm(k = 2, ix = 4, l = 1+1, j = 7)
model_l1.1 <- glm(hf[[2]], family = Gamma(link = "inverse"),
                  data = hf[[1]], na.action = na.omit)
(round((1-model_l1.1[["deviance"]]/model_l1.1[["null.deviance"]])*100,
       digits = 2)) # explained deviance
check_collinearity(model_l1.1) # VIF < 10
summary(model_l1.1)
shapiro.test(model_l1.1$residuals) # p > 0.05 residuals normally distributed
# i.2) Lag 1 - 50th - (river + kd + river:kd)  --------------------------------
hf <- handyf_glm(k = 1, ix = 4, l = 1+1, j = 7)
model_l1.2 <- glm(hf[[2]], family = Gamma(link = "inverse"),
                  data = hf[[1]], na.action = na.omit)
(round((1-model_l1.2[["deviance"]]/model_l1.2[["null.deviance"]])*100,
       digits = 2)) # explained deviance
check_collinearity(model_l1.2) # VIF < 10
summary(model_l1.2)
shapiro.test(model_l1.2$residuals) # p > 0.05 residuals normally distributed
# i.3) Lag 4 - 75th - (enso + kd + sst) ---------------------------------------
hf <- handyf_glm(k = 3, ix = 1, l = 4+1, j = 12)
model_l4.1 <- glm(hf[[2]], family = Gamma(link = "inverse"),
                  data = hf[[1]], na.action = na.omit)
(round((1-model_l4.1[["deviance"]]/model_l4.1[["null.deviance"]])*100,
       digits = 2)) # explained deviance
check_collinearity(model_l4.1) # VIF < 10
summary(model_l4.1)
shapiro.test(model_l4.1$residuals) # p > 0.05 residuals normally distributed
# i.4) Lag 4 - 75th - (enso + kd + wmi) ---------------------------------------
hf <- handyf_glm(k = 3, ix = 1, l = 4+1, j = 13)
model_l4.2 <- glm(hf[[2]], family = Gamma(link = "inverse"),
                  data = hf[[1]], na.action = na.omit)
(round((1-model_l4.2[["deviance"]]/model_l4.2[["null.deviance"]])*100,
       digits = 2)) # explained deviance
check_collinearity(model_l4.2) # VIF < 10
summary(model_l4.2)
shapiro.test(model_l4.2$residuals) # p > 0.05 residuals normally distributed
# i.5) Lag 1 - 75th - (soi + kd) ----------------------------------------------
hf <- handyf_glm(k = 3, ix = 3, l = 1+1, j = 3)
model_l1.3 <- glm(hf[[2]], family = Gamma(link = "inverse"),
                data = hf[[1]], na.action = na.omit)
(round((1-model_l1.3[["deviance"]]/model_l1.3[["null.deviance"]])*100,
       digits = 2)) # explained deviance
check_collinearity(model_l1.3) # VIF < 10
summary(model_l1.3)
shapiro.test(model_l1.3$residuals) # p > 0.05 residuals normally distributed
# i.6) Lag 4 - 75th - (soi + kd + sst) ----------------------------------------
hf <- handyf_glm(k = 3, ix = 3, l = 4+1, j = 12)
model_l4.3 <- glm(hf[[2]], family = Gamma(link = "inverse"),
                  data = hf[[1]], na.action = na.omit)
(round((1-model_l4.3[["deviance"]]/model_l4.3[["null.deviance"]])*100,
              digits = 2)) # explained deviance
check_collinearity(model_l4.3) # VIF < 10
summary(model_l4.3)
shapiro.test(model_l4.3$residuals) # p > 0.05 residuals normally distributed
# i.7) Lag 4 - 75th - (soi + kd + wmi) ----------------------------------------
hf <- handyf_glm(k = 3, ix = 3, l = 4+1, j = 13)
model_l4.4 <- glm(hf[[2]], family = Gamma(link = "inverse"),
                  data = hf[[1]], na.action = na.omit)
(round((1-model_l4.4[["deviance"]]/model_l4.4[["null.deviance"]])*100,
       digits = 2)) # explained deviance
check_collinearity(model_l4.4) # VIF < 10
summary(model_l4.4)
shapiro.test(model_l4.4$residuals) # p > 0.05 residuals normally distributed
# i.8) Lag 9 - 50th - (river + kd + sss + wmmi + river:kd + river:sss) --------
hf <- handyf_glm(k = 2, ix = 4, l = 9+1, j = 28)
model_l9 <- glm(hf[[2]], family = Gamma(link = "inverse"),
               data = hf[[1]], na.action = na.omit)
(round((1-model_l9[["deviance"]]/model_l9[["null.deviance"]])*100,
       digits = 2)) # explained deviance
check_collinearity(model_l9) # VIF < 10
(summary(model_l9)) # river and kd are non-significant
# Stepwise backward prodecure: removing river term
summary(glm(total ~ kd + sss + wmi + river:kd + river:sss,
            family = Gamma(link = "inverse"),
            data = hf[[1]], na.action = na.omit)) # kd is non-significant
# Stepwise backward prodecure: removing kd term
model_stp <- glm(total ~ sss + wmi + river:kd + river:sss,
                  family = Gamma(link = "inverse"),
                  data = hf[[1]], na.action = na.omit)
(round((1-model_stp[["deviance"]]/model_stp[["null.deviance"]])*100,
       digits = 2)) # explained deviance
check_collinearity(model_stp) # VIF < 10
(summary(model_stp)) # all terms are significant
# i.9) Model at lag 10 - 50th - (sam + kd + wmi) ------------------------------
hf <- handyf_glm(k = 2, ix = 2, l = 10+1, j = 13)
hf[[1]] <- na.omit(hf[[1]][1:length(hf[[1]])])
model_l10 <- glm(hf[[2]], family = Gamma(link = "inverse"),
                data = hf[[1]], na.action = na.omit)
DE <- (round((1-model_l10[["deviance"]]/model_l10[["null.deviance"]])*100,
       digits = 2)) # explained deviance
check_collinearity(model_l10) # VIF < 5
(summary(model_l10)) # river and kd are non-significant
# ii) Evaluate the contribution of each variable for GLM selected -------------
m1 <- glm(total ~ kd + wmi, # omit sam 
          family = Gamma(link = "inverse"),
          data = hf[[1]], na.action = na.omit)
m2 <- glm(total ~ sam + wmi, # omit kd
          family = Gamma(link = "inverse"),
          data = hf[[1]], na.action = na.omit)
m3 <- glm(total ~ sam + kd, # omit wmi
          family = Gamma(link = "inverse"),
          data = hf[[1]], na.action = na.omit)
v1de <- de_ev(model_l10, m1)
v2de <- de_ev(model_l10, m2)
v3de <- de_ev(model_l10, m3)
round(v1de/(v1de+v2de+v3de)*DE, digits = 1) # sam
round(v2de/(v1de+v2de+v3de)*DE, digits = 1) # kd
round(v3de/(v1de+v2de+v3de)*DE, digits = 1) # wmi
## Plotting GLM selected ####
model_f <- model_l10
allEffects(model_f)
plot_list <- list()
vrs <-   c("sam", "kd", "wmi") # variables
col_line <- c("#F748A5", "#D55E00", "#000000")
intervals <- list(c(-1.75, 1.75, 0.25), c(-1, 2, 0.25), c(-1.25, 2.25, 0.25))
for(i in 1:3){
  xmin <- intervals[[i]][1]
  xmax <- intervals[[i]][2]
  interval <- intervals[[i]][3]
  p <- plot(glm_effect_single.plot(paste(vrs[i]), model_f,
                                   xmin, xmax, interval,
                                   col_line = col_line[i]))
  plot_list[[i]] <- p + labs(tag = LETTERS[i]) +
    theme(
      plot.tag.position = c(0.14, 0.975),
      plot.tag = element_text(
        face = "bold", size = 20)
    ) 
}
# 
file_name = paste0(file = output , "Fig4_GLM_selected",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
grid.arrange(grobs = plot_list, nrow = 2, ncol = 2)
dev.off()
# 
## Checking RESIDUALS normallity and Plotting RESIDUALS with GRATIA ####
shapiro.test(model_f$residuals) # p > 0.05 residuals are normally distributed
# Observed vs fitted values
obsfit <- appraise(model_f)[[4]]
# Perfect agreement line
model_pal <- lm(data.frame(x = c(10000, 25000), y = c(10000, 25000)))
model_obsfit <- lm(obsfit[["data"]][[2]] ~ obsfit[["data"]][[1]])
correlation <- cor.test(
  (obsfit[["data"]][[2]]), (obsfit[["data"]])[[1]],
  method = "pearson")
r2 <- round(correlation[["estimate"]][["cor"]],2)
pv <- correlation[["p.value"]]
if (pv < 0.001){
  pv <- "< 0.001"
} else {
  pv <- paste("= ", round(pv,3))
}
# Figure
pres <- plot_obs_fit(obsfit, model_pal, 12500, 27500, "A")
file_name = paste0(file = output , "FigS1A_GLM_residuals",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfrow=c(1,1))
pres
dev.off()
# 
# Clear console and environment keeping dataset, packages and functions ####
cat("\014")
ds <- c("path", "data.set", "output", "names_col", "names_lag", "sp")
rm(list = setdiff(ls(), (c(ds, "pk", lsf.str()))))
graphics.off() # Removing figures
gc() # Free unused memory
# ============ GAMs ANALYSES ============================================= ####
## GAMs analyses: TH covariates - TOTAL landings ####
# Diverting the output into a single external file
sink(paste0(file = 'GAMs_Mfurnieri_outputs.txt')) 
# Evaluating different families (fmly) and percentiles (th)
fmly <- list()
fmly[[1]] <- Gamma(link = "inverse")
fmly[[2]] <- inverse.gaussian(link = "1/mu^2")
th <- c('25','50','75')
# i) === Results for each family ===
gams_summary_fm <- list()
for(f in 1:length(fmly)){
  # print(paste0(fmly[[f]][[1]],'(link =',fmly[[f]][[2]], ')'))
  # ii) === Results for each percentile ===
  gams_summary_th <- list()  
  for(k in 1:length(th)){
    # print(th[k])
    # List comprehension to get index (columns number) of variables
    vq <- c(th[k], 'total')
    idx_var <- to_vec(for(i in 1:length(vq))
      str_which(names_col, vq[i], negate = FALSE))
    names_var <- names_col[idx_var] # variable names complete
    # Getting climatic indices by pairs
    ci <- c('soi','sam','enso')
    cixs <- to_list(for(i in 1:length(ci)) 
      c(combn(ci, i, simplify = FALSE)))
    cixss <- list()
    cixss[[1]] <- cixs[[2]] # pairs of climatic indexes
    cixss[[1]][[4]] <- unlist(cixs[[3]]) # all climatic indexes
    names_var_ix <- NA # variable names short for each climatic index
    # iii) === Results for each climatic index and without them ===
    gams_summary_ix <- list() 
    for(ix in 1:length(cixss[[1]])){
      # Getting two of the climatic indices and 'th' landings to omit them
      cith <- paste(paste(cixss[[1]][[ix]], collapse = '|'), 'th', sep = '|')
      names_var_ix = names_var[
        -grep(pattern = cith, names_var)]
      # List comprehension to get index (columns number) of variables
      idx_ss <- to_vec(for(i in 1:length(names_var_ix))
        str_which(names_col, names_var_ix[i], negate = FALSE))
      # Removing tag "th"
      tag <- c(th[k])
      for(i in 1:length(tag)){
        names_var_ix <- gsub(tag[i],'', names_var_ix)
      }
      # print(names_var_ix[1])
      # Covariates and response variable names
      xnames <- head(names_var_ix, -1); ynames <- tail(names_var_ix, 1)
      # GAMs formulas ADDITIVE gams_f <- comb_gam(xnames, ynames, kdf) 
      # GAMs formulas INTERACTION
      inter <- c("river", "kd", "sss") # interaction terms
      gams_f <- comb_gami(xnames, ynames, "4", inter) # getting equations
      if(
        ix < length(cixss[[1]])
      ){
        gams_f <- gams_f[grep(
          pattern = paste(cixss[[1]][[4]], collapse = '|'), (gams_f))]
      } else {
        gams_f <- gams_f
      }
      #  To divert the output into separate external files ordered by family,
      #  percentile and climatic index or river:
      #  - Comment lines 74, 172 and 173
      #  - Uncomment lines 135-138 and 164-165
      # fmly_name <- paste0(fmly[[f]][[1]])
      # sink()
      # sink(paste0(file = path, 'GLMs_Mfurnieri_', fmly_name, '_',
      #             th[k], 'th_', toupper(names_var_ix[1]), '.txt'))
      # iv) === Results for each lag ===
      gams_summary_lg <- list() 
      # colnames(gams_summary_lg) <- c('lag','eq','DE','AICc')
      for(i in 1:length(names_lag)){
        spq <- t(sp[c(idx_ss), , i]) # dataset to each lag
        spq <- as.data.frame(spq)
        # spq$V1 <- as.factor(spq$V1) # climatic index as qualitative variable
        colnames(spq) <- names_var_ix # adding column names to DF
        # v) === Results for each equation ===
        gams_summary_eq <- matrix(NA, ncol = 4, nrow = length(gams_f))
        # output <- matrix(NA, ncol = 4, nrow = 1)
        for(j in 1:length(gams_f)){
          eq <- gams_f[[j]]
          output <- try_gam(eq, spq, i-1, j, fmly[[f]], th[k])
          gams_summary_eq[j,] <- c(
            output$lag, output$eq ,output$DE, output$AICc)
          colnames(gams_summary_eq) <- c('lag','eq','DE','AICc')
          gams_summary_eq <- data.frame(gams_summary_eq)
          gams_summary_eq <- na.omit(gams_summary_eq)
        }
        # Calculate ∆AICc only whether there are valid models
        gams_summary_eq <- unique(gams_summary_eq)
        if (dim(gams_summary_eq)[1] != 0)
        {
          gams_summary_eq$AICc <- (gams_summary_eq$AICc -
                                     min(gams_summary_eq$AICc))
        }
        gams_summary_lg <- rbind(gams_summary_lg, gams_summary_eq)
      }
      # sink() # sets the output back to the console
      # closeAllConnections() # to print on the console again
      gams_summary_ix[[ix]] <- gams_summary_lg
    }
    gams_summary_ix <- setNames(gams_summary_ix, c(rev(ci), 'river'))
    gams_summary_th[[k]] <- gams_summary_ix
  }
  gams_summary_th <- setNames(gams_summary_th, c(th))
  gams_summary_fm[[f]] <- gams_summary_th 
}
sink() # sets the output back to the console
closeAllConnections() # to print on the console again
gams_summary_fm <- setNames(gams_summary_fm, 
                            c(fmly[[1]]$family,fmly[[2]]$family))
save(gams_summary_fm, file = paste0('GAMS_summary_fm.RData'))
# 
## GAMs plot summarized ####
load(paste0("GAMs_summary_fm_FX22.RData"))
# Checking percentiles from output file (.txt) to discard those models with
# high concurvity among variables (worst > 0.5)
ix <- c('enso','sam','soi','river')
# i)  50th --------------------------------------------------------------------
th <- "50"
# === Gamma ===
# Visualize number of GAMs according to each index and their DE values
vis_n_de(gams_summary_fm, fmly[[1]][["family"]], th, ix)
# worst > 0.5 rows 2, 4, 5, 6, 8, 16, 18, 19 have to be dropped
rm_rws <- c(2, 4:6, 8, 16, 18, 19)
lag_v <- vector_de_lag(gams_summary_fm, fmly[[1]][["family"]],
                       th, ix, 'lag')[-c(rm_rws)]
de_v <- vector_de_lag(gams_summary_fm, fmly[[1]][["family"]],
                      th, ix, 'DE')[-c(rm_rws)]
t50g <- table_de_lag(de_v, lag_v, 'gamma') 
# === Inverse gaussian ===
vis_n_de(gams_summary_fm, fmly[[2]][["family"]], th, ix)
# worst > 0.5 rows 1, 2, 4, 12 have to be dropped
rm_rws <- c(1, 2, 4, 12)
lag_v <- vector_de_lag(gams_summary_fm, fmly[[2]][["family"]],
                       th, ix, 'lag')[-c(rm_rws)]
de_v <- vector_de_lag(gams_summary_fm, fmly[[2]][["family"]],
                      th, ix, 'DE')[-c(rm_rws)]
t50i <- table_de_lag(de_v, lag_v, 'inverse.gaussian')
# Merging them into a single DF
cols_num <- c('max','avg','sd')
t50g[cols_num] <- sapply(t50g[cols_num], as.numeric)
t50i[cols_num] <- sapply(t50i[cols_num], as.numeric)
th50 <- rbind (t50g, t50i)
th50$avg <- ifelse(th50$avg == th50$max, NA, th50$avg)
# ii)  25th -------------------------------------------------------------------
th <- "25"
# === Gamma ===
vis_n_de(gams_summary_fm, fmly[[1]][["family"]], th, ix)
# worst > 0.5 rows 1, 2, 3, 4, 6, 20 have to be dropped
rm_rws <- c(1:4, 6, 20)
lag_v <- vector_de_lag(gams_summary_fm, fmly[[1]][["family"]],
                       th, ix, 'lag')[-c(rm_rws)]
de_v <- vector_de_lag(gams_summary_fm, fmly[[1]][["family"]],
                      th, ix, 'DE')[-c(rm_rws)]
t25g <- table_de_lag(de_v, lag_v, 'gamma')  
# === Inverse gaussian ===
vis_n_de(gams_summary_fm, fmly[[2]][["family"]], th, ix)
# worst > 0.5 rows 1, 2, 3, 5 have to be dropped
rm_rws <- c(1:3, 5)
lag_v <-vector_de_lag(gams_summary_fm, fmly[[2]][["family"]],
                      th, ix, 'lag')[-c(rm_rws)]
de_v <- vector_de_lag(gams_summary_fm, fmly[[2]][["family"]],
                      th, ix, 'DE')[-c(rm_rws)]
t25i <- table_de_lag(de_v, lag_v, 'inverse.gaussian')
# Merging them into a single DF
t25g[cols_num] <- sapply(t25g[cols_num], as.numeric)
t25i[cols_num] <- sapply(t25i[cols_num], as.numeric)
th25 <- rbind(t25g, t25i)
th25$avg <- ifelse(th25$avg == th25$max, NA, th25$avg)
# iii)  75th ------------------------------------------------------------------
th <- "75"
# === Gamma ===
# Visualize number of GAMs according to each index and their DE values
vis_n_de(gams_summary_fm, fmly[[1]][["family"]], th, ix)
# worst > 0.5 rows 1, 2, 3, 4, 5, 7, 8, 9, 10,11, 12, 13, 14, 15, 16, 17, 19,
# 20, 21, 23, 35, 37, 38 have to be dropped
rm_rws <- c(1:3, 5, 6:10, 12:14, 16, 26, 27)
lag_v <- vector_de_lag(gams_summary_fm, fmly[[1]][["family"]],
                       th, ix, 'lag')[-c(rm_rws)]
de_v <- vector_de_lag(gams_summary_fm, fmly[[1]][["family"]],
                      th, ix, 'DE')[-c(rm_rws)]
t75g <- table_de_lag(de_v, lag_v, 'gamma')
# === Inverse gaussian ===
vis_n_de(gams_summary_fm, fmly[[2]][["family"]], th, ix)
# worst > 0.5 rows 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
# 28, 30, 31 have to be dropped
rm_rws <- c(1:12, 21, 22)
lag_v <-vector_de_lag(gams_summary_fm, fmly[[2]][["family"]],
                      th, ix, 'lag')[-c(rm_rws)]
de_v <- vector_de_lag(gams_summary_fm, fmly[[2]][["family"]],
                      th, ix, 'DE')[-c(rm_rws)]
t75i <- table_de_lag(de_v, lag_v, 'inverse.gaussian')
# Merging them into a single DF
t75g[cols_num] <- sapply(t75g[cols_num], as.numeric)
t75i[cols_num] <- sapply(t75i[cols_num], as.numeric)
th75 <- rbind (t75g, t75i)
th75$avg <- ifelse(th75$avg == th75$max, NA, th75$avg)
#
# Plot ------------------------------------------------------------------------
p50 <- plot_de_lag(th50) +
  theme(
    legend.direction = "vertical",
    legend.position = c(0.25,0.8),
    legend.key.size = unit(0.75, "cm"),
    legend.text = element_text(size =  12)
    # legend.title = element_text(size = 15, face = "bold")
  ) +
  annotate(geom = "text", x = 2, y = 90, label = "50th",
           color = "black", size = 10, family = "Times") +
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.12, 0.975),
    plot.tag = element_text(
      face = "bold", size = 20)
  ) 
p25 <- plot_de_lag(th25) +
  annotate(geom = "text", x = 2, y = 90, label = "25th",
           color = "black", size = 10, family = "Times") +
  labs(tag = "B") +
  theme(
    plot.tag.position = c(0.12, 0.975),
    plot.tag = element_text(
      face = "bold", size = 20)
  )
p75 <- plot_de_lag(th75) +
  annotate(geom = "text", x = 2, y = 90, label = "75th",
           color = "black", size = 10, family = "Times") +
  labs(tag = "C") +
  theme(
    plot.tag.position = c(0.12, 0.975),
    plot.tag = element_text(
      face = "bold", size = 20)
  )
# Merging plots ---------------------------------------------------------------
file_name = paste0(file = output , "Fig3_GAMs_summarized",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
figure <- ggarrange(p50, p25, p75, ncol = 3, nrow = 1, common.legend = F)
annotate_figure(
  figure, 
  bottom = textGrob(
    expression(paste("Lag (years)")),
    rot = 0, vjust = 0.3,
    gp = gpar(cex = 1.5, fontfamily = "Times" , fontsize = 12)
  ),
  left = textGrob(
    expression(paste("Deviance explained (%)")),
    rot = 90, vjust = 1.3,
    gp = gpar(cex = 1.5, fontfamily = "Times" , fontsize = 12)
  )) 
dev.off()
#
## Checking some models manually ####
# k = 25th: 1 - 50th: 2 - 75th: 3  / ix = ENSO: 1 - SAM: 2 - SOI: 3 - RIVER: 4
# l = lag (from 0 to 10) / j = number of equation (from "gams_summary_fm")
th <- c('25','50','75')
# i.1) Lag 2 - 50th - (sam + sst) ---------------------------------------------
hf <- handyf_gam(k = 2, ix = 2, l = 2+1, j = 5)
model_l2 <- gam(hf[[2]], family = Gamma(link = "inverse"),
                data = hf[[1]], na.action = na.omit, method = "REML")
(round((1-model_l2[["deviance"]]/model_l2[["null.deviance"]])*100,
       digits = 2)) # explained deviance
concurvity(model_l2, full = T) # worst < 0.5
(summary(model_l2)) # sam and sst are significant
shapiro.test(model_l2$residuals) # p > 0.05 residuals are normally distributed
# i.2) Lag 7 - 75th - (sam + kd) ----------------------------------------------
hf <- handyf_gam(k = 3, ix = 2, l = 7+1, j = 3)
model_l7 <- gam(hf[[2]], family = Gamma(link = "inverse"),
                data = hf[[1]], na.action = na.omit, method = "REML")
(round((1-model_l7[["deviance"]]/model_l7[["null.deviance"]])*100,
       digits = 2)) # explained deviance
concurvity(model_l7, full = T) # worst < 0.5
(summary(model_l7)) # sam and kd are significant
shapiro.test(model_l7$residuals) # p > 0.05 residuals are normally distributed
# i.3) Lag 7 - 50th - (sam + kd) ----------------------------------------------
hf <- handyf_gam(k = 2, ix = 2, l = 7+1, j = 3)
model_l7h <- gam(hf[[2]], family = Gamma(link = "inverse"),
                data = hf[[1]], na.action = na.omit, method = "REML")
(round((1-model_l7h[["deviance"]]/model_l7h[["null.deviance"]])*100,
       digits = 2)) # explained deviance
concurvity(model_l7h, full = T) # worst < 0.5
(summary(model_l7h)) # sam and kd are significant
shapiro.test(model_l7h$residuals) # p < 0.05 residuals not normally distributed
# i.4) Lag 9 - 75th - (kd + wmi) ----------------------------------------------
hf <- handyf_gam(k = 3, ix = 4, l = 9+1, j = 13)
hf[[1]] <- na.omit(hf[[1]][1:length(hf[[1]])])
model_l9 <- gam(hf[[2]], family = Gamma(link = "inverse"),
                 data = hf[[1]], na.action = na.omit, method = "REML")
DE <- (round((1-model_l9[["deviance"]]/model_l9[["null.deviance"]])*100,
             digits = 2)) # explained deviance
concurvity(model_l9, full = T) # worst < 0.5
(summary(model_l9)) # kd and wmi are significant
# ii) Evaluate the contribution of each variable for GAM selected -------------
m1 <- gam(total ~ s(kd, k = 4, bs = "cr"), # omit wmi
          family = Gamma(link = "inverse"),
          data = hf[[1]], na.action = na.omit, method = "REML")
m2 <- gam(total ~ s(wmi, k = 4, bs = "cr"), # omit kd
          family = Gamma(link = "inverse"),
          data = hf[[1]], na.action = na.omit, method = "REML")
v1de <- de_ev(model_l9, m1)
v2de <- de_ev(model_l9, m2)
round(v1de/(v1de+v2de)*DE, digits = 1) # kd
round(v2de/(v1de+v2de)*DE, digits = 1) # wmi
# 
## Plotting GAM selected ####
model_f <- model_l9
pf <- plot(model_f)
intercept <- model_f[["coefficients"]][["(Intercept)"]]
cbind(pf[[1]][["x"]], 1/(intercept + pf[[1]][["fit"]]))
cbind(pf[[2]][["x"]], 1/(intercept + pf[[2]][["fit"]]))
plot_list <- list()
vrs <- c("kd", "wmi")
col_line <- c("#D55E00","#000000")
intervals <- list() 
ni <- c(12, 12)
for(i in 1:length(pf)){
  xv <- xaxs(xval = pf, idx = i, precision = 0.05, nintervals = ni[i])
  intervals[i] <- (xv)
}
# Figure
for(i in 1:length(pf)){
  xmin <- intervals[[i]][1]
  xmax <- intervals[[i]][2]
  interval <- intervals[[i]][3]
  p <- plot(gam_effect_single.plot(model_f, intercept, paste(vrs[i]),
                                   xmin, xmax, interval,
                                   col_line = col_line[i]))
  plot_list[[i]] <- p + labs(tag = LETTERS[i]) +
    theme(
      plot.tag.position = c(0.14, 0.975),
      plot.tag = element_text(
        face = "bold", size = 20)
    ) 
}
file_name = paste0(file = output , "Fig5_GAM_selected",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
grid.arrange(grobs = plot_list, nrow = 1, ncol = 2)
dev.off()
## Checking RESIDUALS normallity and Plotting RESIDUALS with GRATIA ####
shapiro.test(model_f$residuals) # p > 0.05 residuals are normally distributed
# Observed vs fitted values
obsfit <- appraise(model_f)[[4]]
# Perfect agreement line
model_pal <- lm(data.frame(x = c(10000, 30000), y = c(10000, 30000)))
model_obsfit <- lm(obsfit[["data"]][[2]] ~ obsfit[["data"]][[1]])
correlation <- cor.test(
  (obsfit[["data"]][[2]]), (obsfit[["data"]])[[1]],
  method = "pearson")
r2 <- round(correlation[["estimate"]][["cor"]],2)
pv <- correlation[["p.value"]]
if (pv < 0.001){
  pv <- "< 0.001"
} else {
  pv <- paste("= ", round(pv,3))
}
# Figure
pres <- plot_obs_fit(obsfit, model_pal, 12500, 27500, 5)
file_name = paste0(file = output , "FigS1B_GAM_residuals",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfrow=c(1,1))
pres
dev.off()
# 