## Script description ####
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Author: Jesus C. Compaire
## Institution: Centro de Investigaciones del Mar y la Atmósfera (CIMA)
## Position: Postdoctoral researcher
## Contact details: jesus.canocompaire@uca.es
## Date created: Nov-2022
##
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## In this script, you will be able to evaluate if the prey abundance ingested 
## by a fish species is related to habitat or particular specimen 
## characteristics by performing a Permutational Multivariate Analysis of
## Variance (PERMANOVA).
## A complete description of PERMANOVA can be checked in Anderson (2001).
## Whilst for a detailed the description of each variable included in the
## dataset check Compaire et al. (2018).
##
## References:
##
## Anderson, M. J. (2001). A new method for non‐parametric multivariate
## analysis of variance. Austral ecology, 26(1), 32-46.
## https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x
##
## Compaire, J. C., Casademont, P., Cabrera, R., Gómez-Cama, C., &
## Soriguer, M. C. (2018). Feeding of Scorpaena porcus (Scorpaenidae) in 
## intertidal rock pools in the Gulf of Cadiz (NE Atlantic). Journal of the
## Marine Biological Association of the United Kingdom, 98(4), 845-853.
## https://doi.org/10.1017/S0025315417000030
##
## -- -- -- -- -- -- -- -- -- -- -- -- --
#
# Clear console and environment ####
cat("\014") 
rm(list=ls())
# Loading package ####
library(vegan)# PERMANOVA is included in the “vegan” package
# Set working directory and load data set from ".CSV" file ####
wd <- setwd("~/Dropbox/GitHub/R/Permanova/")
df <- read.csv("PERMANOVA_dataset.csv", header=TRUE)
attach(df) # To avoid writing "df$column" all the time
str(df) # Displaying internal structure of a data frame in a compact way
# - The first four columns (explanatory variables) are character objects
# - (which are used to represent string values) and the fifth is an integer
# - (contains whole numbers).
# We need to define each explanatory variable as class variable or factor
for (i in 1:length(df)-1){
  df[colnames(df)[i]] <- as.factor(df[ ,i])
  }
# Prey abundance must be defined as an array that contains only numbers. In our
# case, that column is already an object containing only numbers, but if you
# had to modify it, use the next command
df[colnames(length(df))] <- as.matrix(df[ ,length(df)])
# Again, we check the internal structure of the data frame
str(df) 
# - The first four columns are factors and the last one is numeric (integer).
# - Now, we can perform the PERMANOVA test
# PERMANOVA test ####
# Several distance measures are available, in this example we use "Bray-Curtis"
adonis(Prey ~ Depth*Surface, distance = "bray", permutations=9999)
# - The effect of habitat features on prey abundance was not significant for
# - any of the variables (Depth: Pseudo-F = 0.7839, P = 0.537 / Surface: 
# - Pseudo-F = 0.4019, P = 0.722) neither for the interaction between both
# - (Depth x Surface: Pseudo-F = 1.4793; P > 0.2191).
adonis(Prey ~ Sex*SizeClass, distance = "bray", permutations=9999)
# - The effect of specimen characteristics on prey abundance was significant
# - for the interaction between Sex and Size class (Pseudo-F = 2.9659; 
# - P = 0.0208), but not individually (Sex: Pseudo-F = 1.7855, P = 0.143 /
# - SizeClass: Pseudo-F = 0.2713, P = 0.836).
detach(df) 


