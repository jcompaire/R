## Code description ####
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Author: Jesus C. Compaire
## Institution: Centro de Investigaciones del Mar y la Atmósfera (CIMA)
## Position: Postdoctoral researcher
## Contact details: jesus.canocompaire@uca.es
## Date created: Jan-2023
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Code to replicate the statistical analyses and generate figures 
## performed in the manuscript:
## Compaire, J.C., Acha, E.M., Moreira, D. & Simionato C.G. (2023).
## Time series modeling of coastal fishery landings on the Southwestern
## Atlantic shelf
##
## -- -- -- -- -- -- -- -- -- -- -- -- --
#
# Clear console and environment ####
cat("\014") 
rm(list=ls()) 
# Load packages needed for analysis and figures ####
pk <- c("astsa", "car", "comprehenr", "devtools", "FitAR", "forecast",
        "ggplot2", "ggsn", "ggspatial", "ggthemes", "lubridate", "Metrics",
        "ncdf4", "patchwork", "plyr", "dplyr", "raster", "Rcpp", "readxl",
        "rnaturalearth", "sf", "scales", "stats", "tseries", "zoo")
lapply(pk, require, character.only = TRUE)
# Set working directory, invoke functions and load datasets  ####
setwd("~/gdrive/GitHub/R/TimeSeries/") # directory
source("time_series_aucfz_functions.R") # functions
load("time_series_aucfz_data.RData") # data
output <- paste0(getwd(),"/output_figures/") # Directory to put figures
# ============ MAP ======================== ####
aucfz <- as.data.frame(data.set[1]) # Polygon
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
file_name = paste0(file = output , "Fig1_Map",
                   ".pdf", sep="")
pdf(file_name, width=9, height=6, bg = "transparent", colormodel = "cmyk", 
    compress = T)
par(mfcol=c(1,1))
p1 <- map_aucfz(df = bathy, plygn = aucfz, scale = "discrete", clabs)
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
# Clear console and environment keeping dataset, packages and functions ####
cat("\014") 
rm(list = setdiff(ls(), (c("data.set", "output", "pk", lsf.str()))))
gc() # Free unused memory
# ============ CPUE vs COMMERCIAL LANDINGS ======================== ####
df <- as.data.frame(data.set[2])
file_name = paste0(file = output , "Fig2_CPUE_vs_LANDINGS",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
cpue_landings.plot(df)
dev.off()
rm(df)
# ============ BOX PLOTS OF COMMERCIAL LANDINGS =================== ####
df <- list()
for(p in 3:5){
  df[[p-1]] <- data.set[[p]][1:288,] # First 24yr
}
# Combine the list of data frames into one data frame by row
df <- bind_rows(df, .id = "label")
file_name = paste0(file = output , "Fig3_BOXPLOTS_LANDINGS",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
bp_landings.plot(df)
dev.off()
rm(df)
#  
# ============ TIME SERIES ANALYSIS FOR EACH SPECIES ============== ####
# Combine the list of data frames into one data frame by row
df <- bind_rows(data.set[3:5], .id = "label")
# ============ Cynoscion guatucupa - Stripped weakfish ============ ####
#  https://www.fishbase.se/summary/52918
sp <- subset(df, label == 1) # Subset by species
spname <- sp$Scientific_name[1] # Get scientific name
## Visualizing the time series ####
## Saving plot
spname_file <- gsub(". ", "", spname) 
file_name = paste0(file = output,'Fig4a_', spname_file, "_TimeSeries",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
monthly_peaks.plot(sp) +
  plot_annotation(title = spname, subtitle = "A",
                  theme = theme(
                    plot.title = element_text(
                      size = 20, family = "Times", face = "italic",
                      hjust = 0.5, vjust = -7),
                    plot.subtitle = element_text(
                      size = 16, family = "Times", face = "bold",
                      hjust = 0, vjust = -0.25)))
dev.off()
#
## Transform landing data (numeric) to time series ####
x.ts <- ts(sp[c("Landing")],
           start = c(1996,1), end = c(2021,12), frequency = 12)
# 
## We used the first 24 yr (288 months) to fit the model and the last 2 yr
#  (24 months) to evaluate the forecast
x.ts_24yr <- window(x.ts, start = c(1996,1), end = c(2019,12), frequency = 12)
x.ts_2yr <- window(x.ts, start = c(2020,1), end = c(2021,12), frequency = 12)
sp24yr <- sp[1:288,] # Dataframe
# 
## Correlogram of the original series -----------------------------------------
acf2(x.ts_24yr, max.lag = 36, plot = TRUE, main = spname)
#  The peak at 1 year indicates seasonal variation, reflecting a positive linear
#  relationship between pairs of variables (xt, xt+12) separated by 12-months.
#
## Visualizing the time series components (seasonal, trend and irregular) -----
#  Decompose function decomposes a series into the components trend, seasonal
#  effect, and residual
x.ts.decom = decompose(x.ts_24yr, type = "additive", filter = NULL)
## Saving plot
file_name = paste0(file = output,'FigS1a_', spname_file, "_DECOMPLOT",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
dp <- decomp.plot(x = x.ts.decom,
                  main = bquote(paste(
                    "Decomposition of additive time series for ",
                    italic(.(spname)))))
dp + plot_annotation(subtitle = "A",
                     theme = theme(
                       plot.subtitle = element_text(
                         size = 18, family = "Times", face = "bold",
                         hjust = 0.02, vjust = -4)))
dev.off()
#  The plot of the decomposed original time series shows how the trend ranged
#  over time.
#
## Checking homogeneity of variances among years using LEVENE’s test ----------
leveneTest(Landing ~ as.factor(Year), data = sp24yr, center = median)
#  p < .05 variances among years are different -> Transformation is mandatory 
#  to stabilize the variance -> We use the Box-Cox transformation
#
## Transforming the time series -----------------------------------------------
#  First, we calculate lambda,and then we run the transformation
lambda.x <- BoxCox.lambda(x.ts_24yr, method = "guerrero", lower = -2, upper = 2)
x.ts.BoxCox <- tran_boxcox(x.ts_24yr, lambda.x)
#
## Removing trend -------------------------------------------------------------
#  Once the variance was stabilized, we remove the trend differencing the time
#  series to make it stationary. We took differences at lags 1 (non-seasonal)
#  and 12 (seasonal)
diff.x.ts = diff(x.ts.BoxCox, lag = 1) # Differences at lag 1
diff.x.ts.12 = diff(diff.x.ts, lag = 12) # Differences at lag 12
#
#  Augmented Dickey-Fuller test evaluates if the time series has a unit root
tseries::adf.test(diff.x.ts.12, k = 12) # lag order = 12 since data are monthly
#  The Dickey-Fuller unit-root test (p < 0.01) confirms that the time series
#  is stationary
#
## Correlogram of the stationary time series -----------------------------------
#  The correlograms of the autocorrelation function (ACF) and the partial
#  autocorrelation function (PACF) are analyzed on the stationary time series
#  to evaluate the lag orders of the AR and MA processes
## Saving plot
file_name = paste0(file = output, 'Fig5a_', spname_file, "_CORRELOGRAMS",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
correlograms.plot(x = diff.x.ts.12, lags = 48,
                  main = spname,
                  upperindex = "A")
dev.off()
#
## Algorithm to select best ARIMA model ---------------------------------------
#  Fitting different models and selecting the one with the lowest AIC and whose
#  p-value is greater than 0.05 (which means that so residuals are independent)
pi = c(0:3); qi = c(0:2); d = 1; # Non-seasonal ARIMA parameters
Pi = c(0:1); Qi = c(0:1); D = 1; per = 12; # Seasonal ARIMA parameters
run_models(y = x.ts.BoxCox, pi, qi, d, Pi, Qi, D, per)
#
## Model selected -> ARIMA(p,d,q)(P,D,Q)S = ARIMA(1,1,1,0,1,1)12 --------------
pdq = c(1,1,1); PDQ = c(0,1,1); per = 12
best_model <- Arima(x.ts.BoxCox,
                    order = pdq, seasonal = list(order=PDQ, period=per))
## Graphical DIAGNOSIS OF THE RESIDUALS ----------------------------------------
## Saving plot
file_name = paste0(file = output, "Fig6_", spname_file,
                   "_RESIDUAL_DIAGNOSIS", ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
## A) Temporal evolution of residuals ====
res_bm <- best_model$residuals
stdres <- (res_bm)/(sqrt(best_model$sigma2))
pRES <- residuals_ts.plot(x = stdres, upperindex = 'A')
## B) ACF and PACF of residuals ====
pCOR <- correlograms.plot(x = res_bm, lags = 36,
                          main = NULL, upperindex = "B")
## C) Histogram of residuals ====
pval_bm <- Box.test(best_model$residuals,
                    lag = log(length(best_model$residuals)),
                    type = "Ljung-Box")$p.value
pHIS <- residuals_hist.plot(x = stdres, pval = pval_bm, upperindex = 'C')
## D) Scores of the single p-values from the Ljung–Box test ====
boxresult<-as.data.frame(
  LjungBoxTest(res_bm, k = 3, StartLag = 4, lag.max = 36))
pPV <- residuals_pv.plot(boxresult, upperindex = 'D')
#
## Merging each figure into a single plot ====
par(mfcol=c(1,1))
pt <- pRES / pHIS / pCOR + pPV + plot_layout(ncol = 2)
pt +  plot_annotation(title = spname,
                      theme = theme(
                        plot.title = element_text(
                          size = 20, family = "Times",
                          face = "italic",
                          hjust = 0.5, vjust = -0.5)))
dev.off()
## Discussion of DIAGNOSIS OF THE RESIDUALS -----
#  A) Standardized residuals looks like white noise (there is no trend, no
#     outliers and no changing variance across time).
#  B) There is no significant autocorrelation peaks in none of correlograms.
#  C) Residuals are close from normality, and the scores of the p-value from
#     the Ljung–Box test to confirm that they are no significant (p > 0.05).
#  D) p-values are above 0.05 (they are non-significant).
#
## MASE test ------------------------------------------------------------------
#  Evaluating the model performance using the mean absolute scaled error (MASE)
#  Average values of fitted distribution
fit_values <- ts(bxcx(best_model$fitted,lambda.x, 
                      InverseQ = TRUE, type = "BoxCox"),
                 start = c(1996,1), frequency = 12)
MASE <- mase(x.ts_24yr, fit_values, step_size = 1)
#  A MASE value < 1 involves that the actual forecast gives, on average,
#  smaller errors than the one-step errors from a naïve method
# 
## Getting 95% prediction interval for best model ------------------------------
upper <- best_model$fitted + 1.96*sqrt(var(res_bm));
lower <- best_model$fitted - 1.96*sqrt(var(res_bm))
up95fit <- (ts(bxcx(upper,lambda.x, InverseQ = TRUE, type = "BoxCox"),
               start = c(1996,1), frequency = 12))
lo95fit <- (ts(bxcx(lower,lambda.x, InverseQ = TRUE, type = "BoxCox"),
               start = c(1996,1), frequency = 12))
lo95fit[lo95fit < 0] <- 0
# 
## Forecast for 2 yr (24 months) ----------------------------------------------
model_forecast <- forecast(x.ts.BoxCox,
                           h = 48, level = c(95),
                           lambda=NULL, biasadj = NULL, model = best_model)
for_values <- ts(bxcx(
  model_forecast$mean,lambda.x, InverseQ = TRUE, type = "BoxCox"),
  start = c(2020,1), frequency = 12)
##  95 % Prediction interval
up95for <- ts(bxcx(
  model_forecast$upper,lambda.x, InverseQ = TRUE, type = "BoxCox"),
  start = c(2020,1), frequency = 12)
lo95for <- ts(bxcx(
  model_forecast$lower,lambda.x, InverseQ = TRUE, type = "BoxCox"),
  start = c(2020,1), frequency = 12)
lo95for[lo95for < 0] <- 0
#
ff <- cbind(for_values, up95for, lo95for)
write_xlsx(as.data.frame(ff), '~/gdrive/GitHub/R/TimeSeries/dff.xlsx')
## Plotting time series (observed, fitted and forecasted values) --------------
## Saving plot
file_name = paste0(file = output, "Fig9_" ,spname_file, "_FINAL_PLOT",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
line_col <- c('black', "#1B9E77", "darkred") # colour lines
shad_col <- c("#2b4f19", "#C51E1E") # colour shaded area
par(mfcol=c(1,1))
obs_fit_for.plot(x = x.ts, y = fit_values, z = for_values, main = spname,
            line_col, shad_col)
dev.off()
#
## ----------------------------------------------------------------------------
# Clear console and environment keeping dataset, packages and functions ####
cat("\014") 
rm(list = setdiff(ls(), (c("data.set","df", "output", "pk", lsf.str()))))
# ============ Micropogonias furnieri - Whitemouth croaker ============ ####
#  https://www.fishbase.se/summary/7620
sp <- subset(df, label == 2) # Subset by species
spname <- sp$Scientific_name[1] # Get scientific name
## Visualizing the time series ####
## Saving plot
spname_file <- gsub(". ", "", spname) 
file_name = paste0(file = output,'Fig4b_', spname_file, "_TimeSeries",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
monthly_peaks.plot(sp) +
  plot_annotation(title = spname, subtitle = "B",
                  theme = theme(
                    plot.title = element_text(
                      size = 20, family = "Times", face = "italic",
                      hjust = 0.5, vjust = -7),
                    plot.subtitle = element_text(
                      size = 16, family = "Times", face = "bold",
                      hjust = 0, vjust = -0.25)))
dev.off()

## Convert landing data to time series ####
x.ts <- ts(sp[c("Landing")],
           start = c(1996,1), end = c(2021,12), frequency = 12)
# 
## We used the first 24 yr (288 months) to fit the model and the last 2 yr
#  (24 months) to evaluate the forecast
x.ts_24yr <- window(x.ts, start = c(1996,1), end = c(2019,12), frequency = 12)
x.ts_2yr <- window(x.ts, start = c(2020,1), end = c(2021,12), frequency = 12)
x.ts_2yr <- tail(x.ts, n = 24)
sp24yr <- sp[1:288,] # Dataframe
# 
## Correlogram of the original series -----------------------------------------
acf2(x.ts_24yr, max.lag = 36, plot = TRUE, main = spname)
#  The peak at 1 year indicates seasonal variation, reflecting a positive linear
#  relationship between pairs of variables (xt, xt+12) separated by 12-months,
#  and negative for those separated by 6 months.
#
## Visualizing the time series components (seasonal, trend and irregular) -----
#  Decompose function decomposes a series into the components trend, seasonal
#  effect, and residual
x.ts.decom = decompose(x.ts_24yr, type = "additive", filter = NULL)
## Saving plot
file_name = paste0(file = output, "FigS1b_",spname_file, "_DECOMPLOT",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
dp <- decomp.plot(x = x.ts.decom,
                  main = bquote(paste(
                    "Decomposition of additive time series for ",
                    italic(.(spname)))))
dp + plot_annotation(subtitle = "B",
                     theme = theme(
                       plot.subtitle = element_text(
                         size = 18, family = "Times", face = "bold",
                         hjust = 0.02, vjust = -4)))
dev.off()
#  The plot of the decomposed original time series shows how the trend ranged
#  over time.
#
## Checking homogeneity of variances among years using LEVENE’s test ----------
leveneTest(Landing ~ as.factor(Year), data = sp24yr, center = median)
#  p > .05 variances among years are not different -> Transformation is
#  not mandatory 
#
## Removing trend -------------------------------------------------------------
#  Once the variance was stabilized, we remove the trend differencing the time
#  series to make it stationary. We took differences at lags 1 (non-seasonal)
#  and 12 (seasonal)
diff.x.ts = diff(x.ts_24yr, lag = 1) # Differences at lag 1
diff.x.ts.12 = diff(diff.x.ts, lag = 12) # Differences at lag 12
#
#  Augmented Dickey-Fuller test evaluates if the time series has a unit root
tseries::adf.test(diff.x.ts.12, k = 12) # lag order = 12 since data are monthly
#  The Dickey-Fuller unit-root test (p < 0.01) confirms that the time series
#  is stationary
#
## Correlogram of the stationary time series -----------------------------------
#  The correlograms of the autocorrelation function (ACF) and the partial
#  autocorrelation function (PACF) are analyzed on the stationary time series
#  to evaluate the lag orders of the AR and MA processes
## Saving plot
file_name = paste0(file = output, 'Fig5b_', spname_file, "_CORRELOGRAMS",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
correlograms.plot(x = diff.x.ts.12, lags = 48,
                  main = spname, upperindex = "B") # Figure
dev.off()
#
## Algorithm to select best ARIMA model ---------------------------------------
#  Fitting different models and selecting the one with the lowest AIC and whose
#  p-value is greater than 0.05 (which means that so residuals are independent)
pi = c(0:3); qi = c(0:3); d = 1; # Non-seasonal ARIMA parameters
Pi = c(0:4); Qi = c(0:3); D = 1; per = 12; # Seasonal ARIMA parameters
run_models(y = x.ts_24yr, pi, qi, d, Pi, Qi, D, per)
#
## Model selected -> ARIMA(p,d,q)(P,D,Q)S = ARIMA(0,1,3,0,1,1)12 --------------
pdq = c(0,1,3); PDQ = c(0,1,1); per = 12
best_model <- Arima(x.ts_24yr,
                    order = pdq, seasonal = list(order=PDQ, period=per))
## Graphical DIAGNOSIS OF THE RESIDUALS ----------------------------------------
## Saving plot
file_name = paste0(file = output, "Fig7_" ,spname_file,
                   "_RESIDUAL_DIAGNOSIS", ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
## A) Temporal evolution of residuals ====
res_bm <- best_model$residuals
stdres <- (res_bm)/(sqrt(best_model$sigma2))
pRES <- residuals_ts.plot(x = stdres, upperindex = 'A')
## B) ACF and PACF of residuals ====
pCOR <- correlograms.plot(x = res_bm, lags = 36, upperindex = 'B')
## C) Histogram of residuals ====
pval_bm <- Box.test(best_model$residuals,
                    lag = log(length(best_model$residuals)),
                    type = "Ljung-Box")$p.value
pHIS <- residuals_hist.plot(x = stdres, pval = pval_bm, upperindex = 'C')
## D) Scores of the single p-values from the Ljung–Box test ====
boxresult<-as.data.frame(
  LjungBoxTest(res_bm, k = 3, StartLag = 4, lag.max = 36))
pPV <- residuals_pv.plot(boxresult, upperindex = 'D')
#
## Merging each figure into a single plot ====
par(mfcol=c(1,1))
pt <- pRES / pHIS / pCOR + pPV + plot_layout(ncol = 2)
pt +  plot_annotation(title = spname,
                      theme = theme(
                        plot.title = element_text(
                          size = 20, family = "Times",
                          face = "italic",
                          hjust = 0.5, vjust = -0.5)))
dev.off()
## Discussion of DIAGNOSIS OF THE RESIDUALS -----
#  A) Standardized residuals looks like white noise (there is no trend, no
#     outliers and no changing variance across time).
#  B) Less than 5% of spikes outside the confidence interval bounds, 3 peaks
#     in the ACF and 1 in the PACF.
#  C) Residuals are close from normality, and the scores of the p-value from
#     the Ljung–Box test to confirm that they are no significant (p > 0.05).
#  D) p-values are above 0.05 (they are non-significant).
#
## MASE test ------------------------------------------------------------------
#  Evaluating the model performance using the mean absolute scaled error (MASE)
#  Average values of fitted distribution
fit_values <- ts(best_model$fitted,
                 start = c(1996,1), frequency = 12)
MASE <- mase(x.ts_24yr, fit_values, step_size = 1)
#  A MASE value < 1 involves that the actual forecast gives, on average,
#  smaller errors than the one-step errors from a naïve method
# 
## Getting 95% prediction interval for best model ------------------------------
up95fit <- best_model$fitted + 1.96*sqrt(var(res_bm));
lo95fit <- best_model$fitted - 1.96*sqrt(var(res_bm))
lo95fit[lo95fit < 0] <- 0
# 
## Forecast for 2 yr (24 months) ----------------------------------------------
model_forecast <- forecast(x.ts_24yr,
                           h = 48, level = c(95),
                           lambda=NULL, biasadj = NULL, model = best_model)
for_values <- ts(model_forecast$mean, start = c(2020,1), frequency = 12)
##  95 % Prediction interval
up95for <- ts(model_forecast$upper, start = c(2020,1), frequency = 12)
lo95for <- ts(model_forecast$lower, start = c(2020,1), frequency = 12)
lo95for[lo95for < 0] <- 0
#
ff <- cbind(for_values, up95for, lo95for)
write_xlsx(as.data.frame(ff), '~/gdrive/GitHub/R/TimeSeries/dff2.xlsx')
## Plotting time series (observed, fitted and forecasted values) --------------
## Saving plot
file_name = paste0(file = output, "Fig10_",spname_file, "_FINAL_PLOT",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
line_col <- c('black', "#D95F02", "darkred") # colour lines
shad_col <- c("#EA9439", "#C51E1E") # colour shaded area
par(mfcol=c(1,1))
obs_fit_for.plot(x = x.ts, y = fit_values, z = for_values, main = spname,
            line_col, shad_col)
dev.off()
# 
## ----------------------------------------------------------------------------
# Clear console and environment keeping dataset, packages and functions ####
cat("\014") 
rm(list = setdiff(ls(), (c("data.set" ,"df", "output", "pk", lsf.str()))))
# ============ Merluccius hubbsi - Argentine hake ============ ####
#  https://www.fishbase.se/summary/325
sp <- subset(df, label == 3) # Subset by species
spname <- sp$Scientific_name[1] # Get scientific name
## Visualizing the time series ####
## Saving plot
spname_file <- gsub(". ", "", spname) 
file_name = paste0(file = output,'Fig4c_', spname_file, "_TimeSeries",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
monthly_peaks.plot(sp) +
  plot_annotation(title = spname, subtitle = "C",
                  theme = theme(
                    plot.title = element_text(
                      size = 20, family = "Times", face = "italic",
                      hjust = 0.5, vjust = -7),
                    plot.subtitle = element_text(
                      size = 16, family = "Times", face = "bold",
                      hjust = 0, vjust = -0.25)))
dev.off()
## Convert landing data to time series ####
x.ts <- ts(sp[c("Landing")],
           start = c(1996,1), end = c(2021,12), frequency = 12)
# 
## We used the first 24 yr (288 months) to fit the model and the last 2 yr
#  (24 months) to evaluate the forecast
x.ts_24yr <- window(x.ts, start = c(1996,1), end = c(2019,12), frequency = 12)
x.ts_2yr <- window(x.ts, start = c(2020,1), end = c(2021,12), frequency = 12)
sp24yr <- sp[1:288,] # Dataframe
#
## Correlogram of the original series -----------------------------------------
acf2(x.ts_24yr, max.lag = 36, plot = TRUE, main = spname)
#  The peak at 1 year indicates seasonal variation, reflecting a positive linear
#  relationship between pairs of variables (xt, xt+12) separated by 12-months.
#
## Visualizing the time series components (seasonal, trend and irregular) -----
#  Decompose function decomposes a series into the components trend, seasonal
#  effect, and residual
x.ts.decom = decompose(x.ts_24yr, type = "additive", filter = NULL)
## Saving plot
file_name = paste0(file = output, "FigS1c_" ,spname_file, "_DECOMPLOT",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
dp <- decomp.plot(x = x.ts.decom,
                  main = bquote(paste(
                    "Decomposition of additive time series for ",
                    italic(.(spname)))))
dp + plot_annotation(subtitle = "C",
                     theme = theme(
                       plot.subtitle = element_text(
                         size = 18, family = "Times", face = "bold",
                         hjust = 0.02, vjust = -4)))
dev.off()
#  The plot of the decomposed original time series shows a clear negative trend
#  over time.
#
## Checking homogeneity of variances among years using LEVENE’s test ----------
leveneTest(Landing ~ as.factor(Year), data = sp24yr, center = median)
#  p > .05 variances among years are not different -> Transformation is
#  not mandatory 
#
## Removing trend -------------------------------------------------------------
#  Once the variance was stabilized, we remove the trend differencing the time
#  series to make it stationary. We took differences at lags 1 (non-seasonal)
#  and 12 (seasonal)
diff.x.ts = diff(x.ts_24yr, lag = 1) # Differences at lag 1
diff.x.ts.12 = diff(diff.x.ts, lag = 12) # Differences at lag 12
#
#  Augmented Dickey-Fuller test evaluates if the time series has a unit root
tseries::adf.test(diff.x.ts.12, k = 12) # lag order = 12 since data are monthly
#  The Dickey-Fuller unit-root test (p < 0.01) confirms that the time series
#  is stationary
#
## Correlogram of the stationary time series -----------------------------------
#  The correlograms of the autocorrelation function (ACF) and the partial
#  autocorrelation function (PACF) are analyzed on the stationary time series
#  to evaluate the lag orders of the AR and MA processes
## Saving plot
file_name = paste0(file = output, "Fig5c_", spname_file, "_CORRELOGRAMS",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
correlograms.plot(x = diff.x.ts.12, lags = 48,
                  main = spname, upperindex = "C") # Figure
dev.off()
#
## Algorithm to select best ARIMA model ---------------------------------------
#  Fitting different models and selecting the one with the lowest AIC and whose
#  p-value is greater than 0.05 (which means that so residuals are independent)
pi = c(0:4); qi = c(0:1); d = 1; # Non-seasonal ARIMA parameters
Pi = c(0:3); Qi = c(0:1); D = 1; per = 12; # Seasonal ARIMA parameters
run_models(y = x.ts_24yr, pi, qi, d, Pi, Qi, D, per)
#
## Model selected -> ARIMA(p,d,q)(P,D,Q)S = ARIMA(1,1,1,0,1,1)12 --------------
pdq = c(1,1,1); PDQ = c(0,1,1); per = 12
best_model <- Arima(x.ts_24yr,
                    order = pdq, seasonal = list(order=PDQ, period=per))
## Graphical DIAGNOSIS OF THE RESIDUALS ----------------------------------------
## Saving plot
file_name = paste0(file = output, "Fig8_", spname_file,
                   "_RESIDUAL_DIAGNOSIS", ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
## A) Temporal evolution of residuals ====
res_bm <- best_model$residuals
stdres <- (res_bm)/(sqrt(best_model$sigma2))
pRES <- residuals_ts.plot(x = stdres, upperindex = 'A')
## B) ACF and PACF of residuals ====
pCOR <- correlograms.plot(x = res_bm, lags = 36, upperindex = 'B')
## C) Histogram of residuals ====
pval_bm <- Box.test(best_model$residuals,
                    lag = log(length(best_model$residuals)),
                    type = "Ljung-Box")$p.value
pHIS <- residuals_hist.plot(x = stdres, pval = pval_bm, upperindex = 'C')
## D) Scores of the single p-values from the Ljung–Box test ====
boxresult<-as.data.frame(
  LjungBoxTest(res_bm, k = 3, StartLag = 4, lag.max = 36))
pPV <- residuals_pv.plot(boxresult, upperindex = 'D')
#
## Merging each figure into a single plot ====
par(mfcol=c(1,1))
pt <- pRES / pHIS / pCOR + pPV + plot_layout(ncol = 2)
pt +  plot_annotation(title = spname,
                      theme = theme(
                        plot.title = element_text(
                          size = 20, family = "Times",
                          face = "italic",
                          hjust = 0.5, vjust = -0.5)))
dev.off()
## Discussion of DIAGNOSIS OF THE RESIDUALS -----
#  A) Standardized residuals looks like white noise (there is no trend, no
#     outliers and no changing variance across time).
#  B) There is no significant autocorrelation peaks in none of correlograms.
#  C) Residuals are close from normality, and the scores of the p-value from
#     the Ljung–Box test to confirm that they are no significant (p > 0.05).
#  D) p-values are above 0.05 (they are non-significant).
#
## MASE test ------------------------------------------------------------------
#  Evaluating the model performance using the mean absolute scaled error (MASE)
#  Average values of fitted distribution
fit_values <- ts(best_model$fitted,
                 start = c(1996,1), frequency = 12)
MASE <- mase(x.ts_24yr, fit_values, step_size = 1)
#  A MASE value < 1 involves that the actual forecast gives, on average,
#  smaller errors than the one-step errors from a naïve method
# 
## Getting 95% prediction interval for best model ------------------------------
up95fit <- best_model$fitted + 1.96*sqrt(var(res_bm));
lo95fit <- best_model$fitted - 1.96*sqrt(var(res_bm))
lo95fit[lo95fit < 0] <- 0
# 
## Forecast for 2 yr (24 months) ----------------------------------------------
model_forecast <- forecast(x.ts_24yr,
                           h = 48, level = c(95),
                           lambda=NULL, biasadj = NULL, model = best_model)
for_values <- ts(model_forecast$mean, start = c(2020,1), frequency = 12)
##  95 % Prediction interval
up95for <- ts(model_forecast$upper, start = c(2020,1), frequency = 12)
lo95for <- ts(model_forecast$lower, start = c(2020,1), frequency = 12)
lo95for[lo95for < 0] <- 0
#
ff <- cbind(for_values, up95for, lo95for)
write_xlsx(as.data.frame(ff), '~/gdrive/GitHub/R/TimeSeries/dff3.xlsx')
## Plotting time series (observed, fitted and forecasted values) --------------
## Saving plot
file_name = paste0(file = output, "Fig11_", spname_file, "_FINAL_PLOT",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
line_col <- c('black', "#7570B3", "darkred") # colour lines
shad_col <- c(rgb(0,0,1,0.5), "#C51E1E") # colour shaded area
par(mfcol=c(1,1))
obs_fit_for.plot(x = x.ts, y = fit_values, z = for_values, main = spname,
            line_col, shad_col)
dev.off()
# 
## ----------------------------------------------------------------------------
# Argentinian landings ============ ####
rm(df)
df <- bind_rows(data.set[6:8], .id = "label")
## Cynoscion guatucupa ---------------------------------------------------------
sp <- na.omit(subset(df, label == 1)) # Subset by species
spname <- sp$Scientific_name[1]
x.ts <- ts(sp[c("Landing")],
           start = c(1998,1), end = c(2021,12), frequency = 12)
x.ts.decom = decompose(x.ts, type = "additive", filter = NULL)
spname_file <- gsub(". ", "", spname) 
file_name = paste0(file = output, "FigS2a_",spname_file, "_DECOMPLOT",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
dp <- decomp.plot(x = x.ts.decom,
                  main = bquote(paste(
                    "Decomposition of additive time series for ",
                    italic(.(spname)))))
dp + plot_annotation(subtitle = "A",
                     theme = theme(
                       plot.subtitle = element_text(
                         size = 18, family = "Times", face = "bold",
                         hjust = 0.02, vjust = -4)))
dev.off()
#
## Micropogonias furnieri -----------------------------------------------------
sp <- na.omit(subset(df, label == 2)) # Subset by species
spname <- sp$Scientific_name[1]
x.ts <- ts(sp[c("Landing")],
           start = c(1997,1), end = c(2021,12), frequency = 12)
x.ts.decom = decompose(x.ts, type = "additive", filter = NULL)
spname_file <- gsub(". ", "", spname) 
file_name = paste0(file = output, "FigS2b_",spname_file, "_DECOMPLOT",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
dp <- decomp.plot(x = x.ts.decom,
                  main = bquote(paste(
                    "Decomposition of additive time series for ",
                    italic(.(spname)))))
dp + plot_annotation(subtitle = "B",
                     theme = theme(
                       plot.subtitle = element_text(
                         size = 18, family = "Times", face = "bold",
                         hjust = 0.02, vjust = -4)))
dev.off()
#
## Merluccius hubbsi ----------------------------------------------------------
sp <- na.omit(subset(df, label == 3)) # Subset by species
spname <- sp$Scientific_name[1]
x.ts <- ts(sp[c("Landing")],
           start = c(1998,1), end = c(2021,12), frequency = 12)
x.ts.decom = decompose(x.ts, type = "additive", filter = NULL)
spname_file <- gsub(". ", "", spname) 
file_name = paste0(file = output, "FigS2c_",spname_file, "_DECOMPLOT",
                   ".pdf", sep="")
pdf(file_name, width=12, height=6, bg = "transparent", colormodel = "cmyk")
par(mfcol=c(1,1))
dp <- decomp.plot(x = x.ts.decom,
                  main = bquote(paste(
                    "Decomposition of additive time series for ",
                    italic(.(spname)))))
dp + plot_annotation(subtitle = "C",
                     theme = theme(
                       plot.subtitle = element_text(
                         size = 18, family = "Times", face = "bold",
                         hjust = 0.02, vjust = -4)))
dev.off()
#
##
