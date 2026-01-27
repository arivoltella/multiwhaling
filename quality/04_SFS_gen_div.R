#################### 27/01/26 #####################
###### Compute the SFS & genetic diversity ########
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
library(LEA)


#### Read the data #  
VCF1 = readRDS("data/VCF_filtered.RDS")
list_pop <- readRDS("data/list_pop.RDS")


# --------------------------------------------------------------------------------------------------
################ DIVERSITÉ GÉNÉTIQUE et SFS (toutes les pops) ####### ------------------------------
# --------------------------------------------------------------------------------------------------

### Fonctions sources 
source("quality/functions_for_td.r")

data_bin_all <- vcfR2DNAbin(VCF1, extract.indels = TRUE, consensus = FALSE,
                            unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL, verbose = TRUE) 

sfs_tot_all <- as.vector(site.spectrum(data_bin_all, folded = TRUE)) #### sfs folded con funzione PEGAS

#### Plot SFS : 
pdf("plot/SFS/all_pop_sfs.pdf")
barplot(sfs_tot_all) 
dev.off()

# Pour regarder si le SFS dévie du modèle neutre : 
norm_sfs_tot_all <- calcola_normalized_foldedSFS(sfs_tot_all) 

# ADD transformed spectrum of a Kingman pop #
# This is done to check the variation in the SFS class compared to Kingman's expectations
kingman_all <- coal_model(2*length(sfs_tot_all), loci_number = 30000, loci_length = 100, ploidy = 1) + 
  feat_mutation(1) + sumstat_sfs() 
sumstats <- simulate(kingman_all, seed = 23432)
sim_folded_sfs <- fold_dafSFS_coala(sumstats$sfs)
norm_simulated_sfs <- calcola_normalized_foldedSFS(sim_folded_sfs)

######### Plot normalized SFS + modèle neutre : 
pdf("plot/SFS/all_pop_norm_sfs.pdf")
plot(norm_sfs_tot_all, type = "l", ylim = c(0, 0.5))
lines(norm_simulated_sfs, col = 4)
dev.off()

# Métriques de diversité génétique : 
div_all <- unlist(calcola_TD_folded(sfs_tot_all)) 
names(div_all) <- c("TD", "Pi", "W", "S")
div_all
write.csv(div_all, "plot/SFS/all_div.csv")




#--------------------------------------------------------------------------------------------------
########## SFS et DIVERSITÉ GÉNÉTIQUE pop par pop ########## --------------------------------------
#--------------------------------------------------------------------------------------------------


#### Création des listes par population : 
for (i in seq_along(list_pop)) {
  
  # Tri des individus pour les assigner la pop à laquelle on s'intéresse
  pop <- list_pop[i]
  data <- VCF1[,c("FORMAT", unlist(pop))] 
  
  # Calcul du SFS : 
  dati_bin <- vcfR2DNAbin(data, extract.indels = TRUE, consensus = FALSE, 
                          unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL,
                          verbose = TRUE)
  sfs_tot <- as.vector(site.spectrum(dati_bin, folded = TRUE)) 
  
  ######## Plot SFS  : 
  pdf(paste("plot/SFS/", names(list_pop[i]), "_sfs.pdf", sep = ""))
  barplot(sfs_tot) 
  dev.off()
  
  norm_sfs_tot <- calcola_normalized_foldedSFS(sfs_tot) 
  ## transformed SFS to inspect deviation from the neutral model
  ## ADD transformed spectrum of a Kingman pop ########### This is done to check the variation in the SFS class compared to Kingman's expectations
  kingman <- coal_model(2*length(sfs_tot), loci_number = 30000, loci_length = 100, ploidy = 1) + 
    feat_mutation(1) + sumstat_sfs() 
  sumstats <- simulate(kingman, seed = 23432)
  sim_folded_sfs <- fold_dafSFS_coala(sumstats$sfs)
  norm_simulated_sfs <- calcola_normalized_foldedSFS(sim_folded_sfs)
  
  ######### Plot normalized SFS + modèle neutre : 
  pdf(paste("plot/SFS/", names(list_pop[i]), "_norm_sfs.pdf", sep = ""))
  plot(norm_sfs_tot, type = "l", ylim = c(0, 0.5))
  lines(norm_simulated_sfs, col=4)
  dev.off()
  
  ######### Diversité génétique : Tajima's D, thetaPi, thetaS, S
  div_gen <- unlist(calcola_TD_folded(sfs_tot)) 
  names(div_gen) <- c("TD", "Pi", "W", "S")
  write.csv(div_gen, paste("plot/SFS/", names(list_pop[i]), "_div_gen.csv", sep = ""))
}

