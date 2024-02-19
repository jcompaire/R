## Script description ####
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Author: Jesus C. Compaire
## Institution: Centro para el Estudio de Sistemas Marinos (CESIMAR)
## Position: Postdoctoral researcher
## Contact details: jesus.canocompaire@uca.es
## Date created: Jan-2024
##
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## In this script, you will be able to evaluate whether the length and body
## condition, and growth rate of a fish species is influenced by the protection
## status of the shore where inhabits.
## To learn more about this study, you may read the manuscript Compaire et al. 
## (2024). Whilst for a detailed the description of each variable included in
## the dataset check Compaire et al. (2021).
##
## References:
##
## Compaire, J. C, Visintini, N., Soriguer, M. C., Johnson, M. L., Hull, S. &
## Barrett, C. J. (2024). Lipophrys pholis is larger, grows faster and is in
## better condition in protected than in unprotected rocky shores.
## Aquatic Conservation: Marine and Freshwater Ecosystems, 34(2), e4083.
## https://doi.org/10.1002/aqc.4083
##
## Compaire, J. C, Visintini, N., Soriguer, M. C., Johnson, M. L., Hull, S. &
## Barrett, C. J. (2021). Length and weight data of common blenny
## Lipophrys pholis (Blenniidae) caught from tide pools in two contrasting 
## marine provinces in the temperate Northern Atlantic. PANGAEA,
## https://doi.pangaea.de/10.1594/PANGAEA.932955
##
## Note: Small differences between the values of the results presented here and
## those described in Compaire et al., (2024) may arise due to variations in
## the number of decimal places used in the analyses. However, the patterns
## shown here align with those described in the manuscript (such minor 
## differences are indicated with <# ... #>  symbols).
##
## -- -- -- -- -- -- -- -- -- -- -- -- --
#
# Clear console and environment ####
cat("\014") 
rm(list=ls())
# 
# Loading packages ####
pk <- c("dplyr", "emmeans" ,"FSA", "pangaear", "PMCMRplus","readxl", "rstatix")
lapply(pk, require, character.only = TRUE)
# Accessing data hosted in PANGAEA digital repository ####
repo <- pg_data(doi = '10.1594/PANGAEA.932955')
df <- repo[[1]][["data"]]
df <- df[, c(2, 3, 5, 12, 13)] # Subsetting dataframe 
names(df)[c(1, 3, 4, 5)] <- c(
  "Shore", "Date", "TL_cm", "WT_g") # renaming columns
# Adding columns necessary to perform analyses
df$TL_mm <- df$TL_cm*10 # length in mm
df$LogTL <- log10(df$TL_cm) 
df$LogWT <- log10(df$WT_g)
# The body condition is determined using Fulton's condition factor (K)
df$K = (df$WT_g / df$TL_mm^3)*10^5
# Seasons (warm and cold)
df$month_n <- lubridate::month(ymd(df$Date))
df$Season <- ifelse(df$month_n %in% c(6:10), "w", "c")
as.data.frame(df <- subset(df, select = -c(month_n))) # removing columns
# 
# === Table 3 - TL and K ====================================================== 
# Subset by province: "Lusitania" or "Northern European Seas"
subdf <- subset(df, Province == "Lusitania")
# Selecting query variable "TL_cm" or "K"
variable_x <- subdf$TL_cm
factor_x <- as.factor(subdf$Shore) # Query factor convert to character
# Check normality and homoscedasticity
shapiro.test(variable_x) # p < 0.05 we can't assume normality
bartlett.test(variable_x ~ factor_x, data=subdf) # p < 0.05 we can reject
# the null hypothesis that the variance is the same for all treatment groups
# If the result of some of these tests (or both) is p < 0.05 we have to perform
# a non-parametric analysis.
# Kruskal-wallis test to evaluate differences according to type of shore
with(subdf, tapply(variable_x, factor_x, median, na.rm=TRUE))
kruskal.test(variable_x ~ factor_x, data=subdf)
#  Lusitania
#    TL df = 1, H =  1.12, p < 0.291 <# (p = 0.289) #> 
#    K  df = 1, H = 24.17, p < 0.001 <# (H = 24.11) #> 
#  Northern European Seas
#    TL df = 1, H = 77.39, p < 0.001
#    K  df = 1, H = 58.57, p < 0.001 <# (H = 58.11) #>
# Mean, SD, n
with(subdf, tapply(variable_x, factor_x, mean, na.rm=TRUE))
with(subdf, tapply(variable_x, factor_x, sd, na.rm=TRUE))
aggregate(variable_x ~ factor_x, subdf, function(x) sum (x > 0, na.rm = TRUE))
# 
# === Table 4 - Length-Weight Relationships ===================================
# ANCOVA test to assess different growth patterns dividing the population on
# different size classes for fish caught at unprotected shore of NES
subdf <- subset(df, Province == "Northern European Seas" &
                  Shore == "Unprotected")
subdf$size <- with(subdf, ifelse(TL_cm >= 5, 'larger', "smaller"))
#  == Regression line for all population ==
reg_pop <- lm(LogWT ~ LogTL, data=subdf)
summary(reg_pop)
subdf %>%  anova_test(LogWT ~ LogTL * size)
aggregate(TL_cm ~ size, subdf, function(x) sum (x > 0, na.rm = TRUE))
#    There was homogeneity of regression slopes as the interaction term was not
#    statistically significant
#    F(1,147) = 1.353, p = 0.247 <# (F = 1.478, p = 0.226) #>
#    Population should not be separated by size classes.
# 
# ANCOVA test to assess different growth patterns according to shore type for
# each province "Lusitania" or "Northern European Seas"
subdf <- subset(df, Province == "Northern European Seas") 
#  == Regression line for all population ==
reg_pop <- lm(LogWT ~ LogTL, data=subdf)
summary(reg_pop)
subdf %>%  anova_test(LogWT ~ LogTL * Shore)
#  Lusitania
#    There was homogeneity of regression slopes as the interaction term was not
#    statistically significant
#    F(1, 240) = 0.323, p = 0.571 <# F = 0.623, p = 0.431 #>
#    Populations should not be separated by shore.
#  Northern European Seas
#    There was homogeneity of regression slopes as the interaction term was not
#    statistically significant
#    F(1, 199) = 2.494, p = 0.116 <# F = 2.719, p = 0.101 #>
#    Populations should not be separated by shore.
# Despite there is no difference, we are going to calculate the parameters
# according to type of shore # "Protected" or "Unprotected"
# == Regresion line according to shore type ==
subdf <- subset(subdf, Shore == "Unprotected")
# Minimum and maximum values of length and weight 
round((c(min(subdf$TL_cm), max(subdf$TL_cm))),1)
round((c(min(subdf$WT_g), max(subdf$WT_g))),2)
# Regression
reg_sho <- lm(LogWT ~ LogTL, data=subdf)
summary(reg_sho)
reg_sho[["coefficients"]][["LogTL"]] #  "b" parameter
# Access to intercept value to obtain inverse log of "a"
name_list <- names(reg_sho); # list the names of model variables
intercept_v <- coefficients(reg_sho)["(Intercept)"]; # select intercept value
10^intercept_v; # Inverse log value of 'a'
# Confidence Interval
CI <- confint(reg_sho)
10^CI[1,] # Inverse log values to 95% CI "a"
CI[2,] # Values to 95% CI "b"
#  Lusitania - Protected
#    b = 3.166 (3.097 - 3.235) <# b = 3.173 (3.104 - 3.241) #>
#    a = .0081 (.0073 - .0090) <# a = .0080 (.0072 - .0089) #>
#  Lusitania - Unprotected
#    b = 3.131 (3.032 - 3.231) <# b = 3.125 (3.025 - 3.224) #>
#    a = .0078 (.0067 - .0090) <# a = .0079 (.0068 - .0091) #>
#  Northern European Seas - Protected
#    b = 3.109 (2.978 - 3.241) <# b = 3.110 (2.978 - 3.241) #>
#    a = .0101 (.0077 - .0134) <# a = .0101 (.0077 - .0133) #>
#  Northern European Seas - Unprotected
#    b = 2.954 (2.874 - 3.034) <# b = 2.949 (2.870 - 3.029) #>
#    a = .0108 (.0096 - .0121) <# a = .0109 (.0097 - .0122) #>
# 
# ANCOVA test to assess different growth patterns according to the season for
# each shore and province
subdf <- subset(df, Province == "Northern European Seas" &
                  Shore == "Protected")
subdf %>%  anova_test(LogWT ~ LogTL * Season)
#  Lusitania - Protected
#    There was homogeneity of regression slopes as the interaction term was not
#    statistically significant
#    F(1, 97) = 0.914, p = 0.341 <# F = 1.603, p = 0.209 #>
#    Population should not be separated by season.
#  Lusitania - Unprotected
#    There was homogeneity of regression slopes as the interaction term was not
#    statistically significant
#    F(1, 139) = 1.303, p = 0.256 <# F = 1.912, p = 0.169 #>
#    Population should not be separated by season.
#  Northern European Seas - Protected
#    There was homogeneity of regression slopes as the interaction term was not
#    statistically significant
#    F(1, 48) = 1.912, p = 0.173 <# F = 1.908, p = 0.174 #>
#    Population should not be separated by season.
#  Northern European Seas - Unprotected
#    There was homogeneity of regression slopes as the interaction term was not
#    statistically significant
#    F(1, 147) = 2.215, p = 0.139 <# F = 2.070, p = 0.152 #>
#    Population should not be separated by season.
# 
