#Cargar el fichero

Prey<-read.csv("/Users/Worm/Dropbox/RNM243/Personales/Mila_Gsus/FuturosPapers/Paper_Sporcus_Dieta/4_JMBA/PERMANOVA/PERMANOVA_Prey.csv", header=TRUE)

attach(Prey)

Prey

#Romper el fichero y generar las variables

Depth=as.factor(Prey[ ,1])

Surface=as.factor(Prey[ ,2])

Physiography=as.factor(Prey[ ,3])

Season=as.factor(Prey[ ,4])

SizeClass=as.factor(Prey[ ,5])

Sex=as.factor(Prey[ ,6])

Prey=as.matrix(Prey[ ,7])

#Instalar el paquete "vegan" para hacer el Permutational Multivariate Analysis of Variance.
#You have to choose a distance measure. For regular measurements, that would usually be "euclidean", but for species community data, "bray" (the Bray-Curtis distance) is more appropriate. The following will give you an output exactly equivalent to an ANOVA table:

library(vegan)

adonis(Prey~Depth, distance="bray")

adonis(Prey~Surface, distance="bray")

adonis(Prey~Physiography, distance="bray")

adonis(Prey~Season, distance="bray")

adonis(Prey~SizeClass, distance="bray")

adonis(Prey~Sex, distance="bray")

adonis(Prey~Depth*Surface, distance="bray",permutations=9999)

adonis(Prey~Sex*SizeClass, distance="bray",permutations=9999)


