#################### 24/03/26 #####################
############ Sliding window SFS : #################
############# Achille RIVOLTELLA ##################

# Libraries : 
library(vcfR)
library(tidyverse)
library(rlist)
library(pegas)
library(ggplot2)
library(adegenet)
library(patchwork)

#######################################################################################################
# ARGUMENTS à mettre en entrée :                                                                      #
#         - VCF avec SNPs ; vecteur avec toutes les positions ;                                       #
#         - objet (liste) contenant les noms des pops et les individus associés                       #
#         - n° chromosome                                                                             #
#         - taille de fenêtre et pas                                                                  #  
#                                                                                                     #
# SORTIES de la fonction :                                                                            #
#         - SFS par pop et par fenêtre (sfs_spectre_final)                                            #
#         - Dataset des indices de diversité génétique par pop et par fenêtre (div_gen_final)         #
#######################################################################################################

#### ARGUMENTS :  
args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]
#chrom <- 21

window <- args[2]
slide <- args[3]
#window <- 500000
#slide <- 100000

# Data : 
VCF <- read.vcfR(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/VCF_filtered_", 
                       chrom, ".vcf.gz", sep = ""))

VectPosFiltered <- readRDS(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/vec_pos/vec_pos_all_", 
                                  chrom, ".RDS", sep = ""))                 # Toutes les positions du chromosomes 

list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/data/list_pop.RDS")

# ATTENTION : IMPORTER LES FONCTIONS NÉCESSAIRES 
source("/shared/projects/multiwhaling/multiwhaling/whole_genome/fonctions/fonctions.R")

#### FILTRER LES NA POUR ANALYSE SFS ####

################################ SLIDING WINDOW ################################

fst_sliding_window <- function(VCF, list_pop, window, slide){
  
  # Préparation des bornes à chaque itération : ---------------------------------------
  vec_pos <- getPOS(VCF)
  low_bound <- seq(1, max(vec_pos) + slide, slide)
  upper_bound <- seq(window, max(vec_pos) + slide + window, slide)
  positions_plot <- (upper_bound + low_bound)/2
  n_pop <- length(list_pop)
  
  
  #### Pour chaque fenêtre : ----------------------------------------------------------
  sfs_spectre_sliding <- c()
  fenetre_all <- c()
  
  for (i in 1:length(low_bound)) {
    fenetre <- which(vec_pos >= low_bound[i] & vec_pos <= upper_bound[i])
    fenetre_all[i] <- length(which(VectPosFiltered >= low_bound[i] & VectPosFiltered <= upper_bound[i])) # Nombre total de positions après filtres par fenêtres
    
    # ---- #### ... s'il y a des SNPs dans la fenêtre de positions ...  ---------------
    if (length(fenetre) > 0) {
      data_fenetre <- VCF1[fenetre,]         # On sélectionne cette fenêtre de SNPs
      
      # ---- #### Sélectionner 1 population ------------ 
      temp_spectre <- c()
      for (z in 1:length(list_pop)){
        
        data_temp <- data_fenetre[, c("FORMAT",list_pop[[z]])]
        
        # ------- #### Filtrer les NA ... -----------
        genotypes <- extract.gt(data_temp, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
        NAs <- rowSums(genotypes == "./.")
        data_temp <- subset(data_temp, NAs < 1)
        
        # ------- #### Et le SFS folded ... --------- # Sur le VCF filtré par : fenêtre, population et les NA 
        data_bin <- vcfR2DNAbin(data_temp, extract.indels = TRUE, 
                                consensus = FALSE,  unphased_as_NA = FALSE, ref.seq = NULL, 
                                start.pos = NULL, verbose = TRUE) 
        sfs_tot <- site.spectrum(data_bin, folded = TRUE)
        temp_spectre <- c(temp_spectre, sfs_tot)
      }
      sfs_spectre_sliding <- rbind(sfs_spectre_sliding,temp_spectre) 
      # À chaque itération, ajoute le SFS de toutes les pops côtes à côtes et répète cela pour chaque fenêtre
    }
    # ------- #### Sinon on met des NaN : ------------------------------------------
    else {
      sfs_spectre_sliding <- rbind(sfs_spectre_sliding, rep(NaN,length(unlist(list_pop))))
    }
  }
  
  # ------------------------------------------------------------------------------
  #################### Calculs des indices de diversité ##########################
  # ------------------------------------------------------------------------------
  
  # Mise en forme des données : 
  positions_plot <- (upper_bound + low_bound)/2
  sfs_spectre_final <- as.data.frame(na.omit(cbind(positions_plot, sfs_spectre_sliding)))
  write_csv(sfs_spectre_final, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/SFS/", chrom, "_sfs_spectre_sliding_", window, "_", slide, ".csv", sep = ""))
  
  # Calculer SFS norm par pop : 
  a <- 0
  div_gen_final <- c()
  
  for (i in seq_along(list_pop)) {
    
    # ----------  #### Sélection de la pop : 
    n_cols <- length(list_pop[[i]])      # Nb d'individus appartenant à la pop i
    col_range <- (a + 1):(a + n_cols)    # Sélectionne les individus (colonnes) à partir de la colonne suivant la dernière de la pop précédente ... dernier de la pop
    
    # ----------  #### Calcul des indices de diversité : TD, Pi, W et S
    matrice_sfs_bypop <- sfs_spectre_sliding[, col_range]
    qq = apply(matrice_sfs_bypop, 1,calcola_TD_folded)          # Calcul du SFS ligne par ligne          
    divgen = matrix(unlist(qq), ncol= 4, nrow = nrow(matrice_sfs_bypop), byrow = T)
    divgen_spectre <- as_tibble(cbind(positions_plot, divgen), .name_repair = "unique")
    colnames(divgen_spectre) <- c("positions", "TD", "Pi", "W", "S")
    
    
    # ---------- #### Faire le scaling par le nombre de positions comprises dans chaque fenêtre 
    divgen_spectre <- divgen_spectre |>
      mutate(Pi = Pi/fenetre_all,            # Pour les deux estimateurs, on scale par le nb total de positions séquencées (SNP & NPP)
             W = W/fenetre_all,              
             Pop = paste(names(list_pop[i])))              
    
    div_gen_final <- rbind(div_gen_final, divgen_spectre)

    # ---------- #### On passe à la pop suivante : 
    a <- a + n_cols
  }
  saveRDS(div_gen_final, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/SFS/", chrom, "_div_gen_final_", window, "_", slide, ".csv", sep = ""))
}

