#################### 19/03/26 #####################
###### Compute the SFS & genetic diversity ########
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(ggplot2)
library(tidyverse)
library(pegas)
library(coala)

#######################################################################################################
# ARGUMENTS à mettre en entrée :                                                                      #
#         - VCF avec SNPs ; VCF avec toutes les positions ;                                           #
#         - objet (liste) contenant les noms des pops et les individus associés                       #
#         - n° chromosome                                                                             #
#                                                                                                     #  
#                                                                                                     #
# SORTIES de la fonction :                                                                            #
#         - Objets contenant SFS et SFS norm par chrom et par pop                                     #
#         - plot de SFS par chromosome et par population                                              #
#         - plot de SFS normalisé + modèle neutre par chrom et par pop                                #
#         - métriques de diversité génétique (TD, Pi, W, S) par chrom et par pop                      #
#######################################################################################################


#### ARGUMENTs : ---------------------------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]
#chrom <- 21

VCF_all = read.vcfR(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/VCF_all_filtered_", 
                           chrom, ".vcf.gz", sep = ""))  
VCF = read.vcfR(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/VCF_filtered_", 
                        chrom, ".vcf.gz", sep = ""))         
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/list_pop.RDS")

### Fonctions sources /!\ À MODIFIER :
source("/shared/projects/multiwhaling/multiwhaling/whole_genome/fonctions/fonctions.R")


SFS_diversity <- function(VCF1, VCF_all, list_pop){     # -----------------------------------------# 
  
  #list_pop$all <- unlist(list_pop, use.names = F)
  
  sfs_pop <- c()
  norm_sfs <- c()
  
  # Sauvegarde du vecteur de position total : (pour sliding window SFS)
  vec_pos_all <- getPOS(VCF_all)
  saveRDS(vec_pos_all, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/vec_pos/vec_pos_all_", 
                             chrom, ".RDS", sep = ""))

  #-------------------------------------------------------------------------------------------------#
  ############################# SFS et DIVERSITÉ GÉNÉTIQUE pop par pop ##############################
  #-------------------------------------------------------------------------------------------------#
  
  #### Création des SFS population par population : 
  for (i in seq_along(list_pop)) {
    
    ####  Tri par population : ----------------------------------------------------------------
    pop <- list_pop[i]
    data <- VCF1[,c("FORMAT", unlist(pop))] 
    data_all <- VCF_all[,c("FORMAT", unlist(pop))]      # On garde les individus d'une pop
    
    ####  Filtre des NAs : --------------------------------------------------------------------
    genotypes <- extract.gt(data, element = "GT", mask = FALSE, as.numeric=F,
                            return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, 
                            convertNA = FALSE)
    NAs <- rowSums(genotypes == "./.")
    data <- subset(data, NAs < 1)
    
    genotypes <- extract.gt(data_all, element = "GT", mask = FALSE, as.numeric=F,
                            return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, 
                            convertNA = FALSE)
    NAs <- rowSums(genotypes == "./.")
    data_all <- subset(data_all, NAs < 1)
    positions_all <- getPOS(data_all)
    
    #saveRDS(positions_all, paste("/shared/projects/multiwhaling/multiwhaling/data/vect_position/pos_", names(list_pop[i]), ".RDS", sep = ""))
    
    #### Calcul du SFS : ----------------------------------------------------------------------
    dati_bin <- vcfR2DNAbin(data, extract.indels = TRUE, consensus = FALSE, 
                            unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL,
                            verbose = TRUE)
    sfs_tot <- as.vector(site.spectrum(dati_bin, folded = TRUE)) 
    sfs_pop <- c(sfs_pop, list(sfs_tot))
    
    #### Plot SFS : 
    pdf(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/04_SFS/", names(list_pop[i]), "_SFS_", chrom, ".pdf", sep = ""))
    barplot(sfs_tot, legend.text = paste(names(list_pop[i]))) 
    dev.off()
    
    #### Normalized SFS + modèle neutre : -----------------------------------------------------
    norm_sfs_tot <- calcola_normalized_foldedSFS(sfs_tot) 
    norm_sfs <- c(norm_sfs, list(norm_sfs_tot))
    
    ## transformed SFS to inspect deviation from the neutral model
    ## ADD transformed spectrum of a Kingman pop ########### This is done to check the variation in the SFS class compared to Kingman's expectations
    kingman <- coal_model(2*length(sfs_tot), loci_number = 30000, loci_length = 100, ploidy = 1) + 
      feat_mutation(1) + sumstat_sfs() 
    sumstats <- simulate(kingman, seed = 23432)
    sim_folded_sfs <- fold_dafSFS_coala(sumstats$sfs)
    norm_simulated_sfs <- calcola_normalized_foldedSFS(sim_folded_sfs)
    
    #### Plot normalized SFS + modèle neutre : 
    pdf(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/04_SFS/", names(list_pop[i]), "_SFS_norm_", chrom, ".pdf", sep = ""))
    plot(norm_sfs_tot, type = "l", ylim = c(0, 0.5))
    lines(norm_simulated_sfs, col=4)
    title(main = paste(names(list_pop[i])))
    dev.off()
    
    #### Diversité génétique : Tajima's D, thetaPi, thetaS, S : -------------------------------
    div_gen <- unlist(calcola_TD_folded(sfs_tot)) 
    names(div_gen) <- c("TD", "Pi", "W", "S")
    div_gen <- as_tibble(t(as.data.frame(div_gen))) |>
      mutate(Pi = Pi/length(positions_all),                       # Pour les deux estimateurs, on scale par le nb total 
             W = W/length(positions_all))                         # de positions séquencées (SNP & NPP)
    write.csv(div_gen, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/04_SFS/", names(list_pop[i]), "_div_gen_wh_", chrom, ".csv", sep = ""))
  }
  saveRDS(sfs_pop, "/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/04_SFS/SFS_all_pop.RDS")
  saveRDS(norm_sfs, "/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/04_SFS/SFS_norm_all_pop.RDS")
}
##################################################################################################################


SFS_diversity(VCF, VCF_all, list_pop)



