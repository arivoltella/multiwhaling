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
VCF_all = readRDS("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_all.RDS")  # On utilise le VCF avec toutes les positions pour scaler les indices de div génétique
VCF1 = readRDS("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered.RDS")         # VCF filtré auquel on va enlever les NA restants par pop
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/list_pop2.RDS")

list_pop$all <- unlist(list_pop, use.names = F)

### Fonctions sources 
source("/shared/projects/multiwhaling/multiwhaling/quality/functions_for_td.r")

#--------------------------------------------------------------------------------------------------
########## SFS et DIVERSITÉ GÉNÉTIQUE pop par pop ########## --------------------------------------
#--------------------------------------------------------------------------------------------------

#### Création des SFS population par population : 
for (i in seq_along(list_pop)) {
  
  ########  Tri par population : 
  pop <- list_pop[i]
  data <- VCF1[,c("FORMAT", unlist(pop))] 
  data_all <- VCF_all[,c("FORMAT", unlist(pop))]        # On garde les individus d'une pop
  
  ########  Filtre des NAs : 
  genotypes <- extract.gt(data, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  NAs <- rowSums(genotypes == "./.")
  data <- subset(data, NAs < 1)
  
  genotypes <- extract.gt(data_all, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  NAs <- rowSums(genotypes == "./.")
  data_all <- subset(data_all, NAs < 1)
  positions_all <- getPOS(data_all)
  
  saveRDS(positions_all, paste("/shared/projects/multiwhaling/multiwhaling/data/vect_position/pos_", names(list_pop[i]), ".RDS", sep = ""))
  
  ######## Calcul du SFS : 
  dati_bin <- vcfR2DNAbin(data, extract.indels = TRUE, consensus = FALSE, 
                          unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL,
                          verbose = TRUE)
  sfs_tot <- as.vector(site.spectrum(dati_bin, folded = TRUE)) 
  
  ######## Plot SFS : 
  pdf(paste("/shared/projects/multiwhaling/multiwhaling/plot/SFS/", names(list_pop[i]), "_sfs.pdf", sep = ""))
  barplot(sfs_tot, legend.text = paste(names(list_pop[i]))) 
  dev.off()
  
  ######## Normalized SFS + modèle neutre :
  norm_sfs_tot <- calcola_normalized_foldedSFS(sfs_tot) 
  ## transformed SFS to inspect deviation from the neutral model
  ## ADD transformed spectrum of a Kingman pop ########### This is done to check the variation in the SFS class compared to Kingman's expectations
  kingman <- coal_model(2*length(sfs_tot), loci_number = 30000, loci_length = 100, ploidy = 1) + 
    feat_mutation(1) + sumstat_sfs() 
  sumstats <- simulate(kingman, seed = 23432)
  sim_folded_sfs <- fold_dafSFS_coala(sumstats$sfs)
  norm_simulated_sfs <- calcola_normalized_foldedSFS(sim_folded_sfs)
  
  ######## Plot normalized SFS + modèle neutre : 
  pdf(paste("/shared/projects/multiwhaling/multiwhaling/plot/SFS/", names(list_pop[i]), "_norm_sfs.pdf", sep = ""))
  plot(norm_sfs_tot, type = "l", ylim = c(0, 0.5))
  lines(norm_simulated_sfs, col=4)
  title(main = paste(names(list_pop[i])))
  dev.off()
  
  ######### Diversité génétique : Tajima's D, thetaPi, thetaS, S
  div_gen <- unlist(calcola_TD_folded(sfs_tot)) 
  names(div_gen) <- c("TD", "Pi", "W", "S")
  div_gen <- as_tibble(t(as.data.frame(div_gen))) |>
    mutate(Pi = Pi/length(positions_all),      # Pour les deux estimateurs, on scale par le nb total 
           W = W/length(positions_all))        # de positions séquencées (SNP & NPP)
  write.csv(div_gen, paste("/shared/projects/multiwhaling/multiwhaling/plot/SFS/", names(list_pop[i]), "_div_gen.csv", sep = ""))
}

