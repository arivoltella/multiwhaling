################### 27/01/26 ######################
########## sNMF clustering analysis ###############
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)
library(pegas)
library(coala)
library(MASS)


#### Read the data #  
VCF1 = readRDS("data/VCF_filtered.RDS")
list_pop <- readRDS("data/list_pop.RDS")

names_ind <- enframe(list_pop, name = "Population", value = "Individu") |>
  unnest(Individu)

# ------------------------------------------------------------------------------------------
################## sNMF clustering ############## ------------------------------------------
# ------------------------------------------------------------------------------------------

whale_snmf <- snmf(input.file= "data/VCF_filtered.geno", K= 1:6, repetitions = 3, project="new", entropy=T)
#show(whale_snmf)
summary(whale_snmf)

pdf(file = "plot/Clustering/entropy_K.pdf")
plot(whale_snmf, col = "blue", pch = 19, cex = 1.2) ### Identifier le K avec l'entropie la plus faible
dev.off()

best <- which.min(cross.entropy(whale_snmf, K = 6))
my.colors <- c("royalblue", "seagreen4", "peru", "khaki", "indianred2", "violet")

adm_coeff <- Q(whale_snmf, K = 4, run = best) 

pdf(file = "plot/Clustering/structure_plot.pdf")
barplot(t(adm_coeff), names = names_ind$Individu, las = 2, cex.names = 0.6, border = NA, space = 0, 
        col = my.colors, 
        xlab = "Individuals", 
        ylab = "Ancestry proportions", 
        main = "Ancestry matrix")
dev.off()
