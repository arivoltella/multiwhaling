################### 12/03/26 ######################
############### scan PCA analysis #################
############## Achille RIVOLTELLA #################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)
library(adegenet)
library(hierfstat)

#######################################################################################################
# ARGUMENTS à mettre en entrée :                                                                      #
#         - VCF avec SNPs ; vecteur toutes les positions (SNP & NPP) ;                                #
#         - objets (liste et tibble) contenant les noms des pops et les individus associés            #
#         - n° chromosome                                                                             #
#         - Tailles de fenêtre et pas à modifier directement dans le script                           #
#                                                                                                     #
# SORTIES de la fonction :                                                                            #
#         - VCF transformé en .geno (par chromosome)                                                  #
#         - Scan d'ACP sur les SNPs : dataset et plot variance PC1 & PC2 le long du chr               #
#         - Scan d'ACP sur les positions :                                                            #
#                             - dataset et plot variance PC1 & 2 le long du chr                       #
#                             - dataset et plot coord PC1 & 2 des individus le long du chr            # 
#######################################################################################################

#### ARGUMENTS :  
args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]
#chrom <- 21

VCF = read.vcfR(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/VCF_filtered_", 
                      chrom, ".vcf.gz", sep = ""))   
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/list_pop.RDS")
list_pop <- list_pop[1:6]
names_ind <- readRDS("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/names_ind.RDS")
VectPosFiltered <- readRDS(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/vec_pos_all_", 
                                 chrom, ".RDS", sep = ""))         


PCA_scan <- function(VCF1, vec_pos, names_ind){
  
  #### Initialization : objets re-used 
  vectNom <- c(names(VCF1@gt[1,])[-1]) 
  vec_snps_pos <- getPOS(VCF1)
  
  #------------------------------------------------------------------------------------------------------
  ################################# PCA sliding window scan on the SNPs ################################# 
  ############################ Representing PC1 & PC2 variance for each window ##########################
  #------------------------------------------------------------------------------------------------------
  
  # Sliding window size : 
  n_snp <- 300
  slide <- 100
  
  snp_total <- length(vec_snps_pos)
  
  low_bound <- seq(1, snp_total, slide)
  upper_bound <- low_bound + n_snp
  
  # Initialize : 
  PC1 <- c()
  PC2 <- c()
  low_f <- c()
  upper_f <- c()
  
  #### For each window : 
  for (i in 1:length(low_bound)) {              
    if (upper_bound[i] <= snp_total) {         # If not yet the end of the VCF 
      
      # ---- #### Convert the VCF to genind ...
      data_moment <- VCF1[low_bound[i]:upper_bound[i],]
      genotype <- extract.gt(data_moment, element = "GT", mask = FALSE, as.numeric=F,
                             return.alleles = FALSE, IDtoRowNames = TRUE, 
                             extract = TRUE, convertNA = FALSE)
      
      geno_temp <- t(genotype)                 #### On transpose car objet genind a besoin d'une matrice 
      #### avec individus en lignes et SNPs en colonnes
      for (j in 1:ncol(geno_temp)){           #### ./. --> NA
        geno_temp[geno_temp[,j]=="./.",j] <- NA 
      }
      mat_data <- df2genind(geno_temp, ploidy = 2, sep="/", pop = vectNom) #### creo oggetto genind dove ho specificato anche le pop di appartenenza
      
      # ---- #### Compute the PCA ... 
      whale_pca <- indpca(mat_data)
      new_eigen <- whale_pca$ipca$eig / sum(whale_pca$ipca$eig)
      PC1 <- c(PC1,new_eigen[1])
      PC2 <- c(PC2,new_eigen[2])      # Récupérer les % de variance des 2 premières PC
      
      low_f <- c(low_f,vec_snps_pos[low_bound[i]])
      upper_f <- c(upper_f,vec_snps_pos[upper_bound[i]])
    }
  }
  
  finale <- cbind(low_f, upper_f, PC1, PC2)  
  
  finale <- as.data.frame(finale) |>
    mutate(position = (low_f + upper_f)/2)
  saveRDS(finale, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/PCA/", chrom, "_scan_pca_snps.RDS", sep = ""))
  
  pdf(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/PCA/", chrom, "_scan_pca_snps.pdf", sep = ""))
  finale |>
    ggplot() + 
    geom_point(aes(x = position, y = PC1)) +
    geom_point(aes(x = position, y = PC2), color = 'red') +
    scale_y_continuous(limits = c(0, 0.2))
  dev.off()
  

  #------------------------------------------------------------------------------------------------------
  ################################ PCA sliding window scan on the positions #############################
  ############################# Representing PC1 & PC2 variance for each window #########################
  ################################### and PC1 & PC2 coords for each indiv ###############################
  #------------------------------------------------------------------------------------------------------
  
  # Sliding window size : ---------------------------
  window <- 500000
  slide <- 100000
  
  vec_snps_pos <- getPOS(VCF1)
  
  low_bound <- seq(1, max(VectPosFiltered) + slide, slide)
  upper_bound <- seq(window, max(VectPosFiltered) + window + slide, slide)
  
  snp_total <- last(VectPosFiltered)
  
  # Initialize : 
  final_acp_var1_2 <- c()
  final_acp_coord1_2 <- c()

  
  #### For each window : ---------------------------------------------------------------------
  for (i in 1:length(low_bound)) {              
    if (upper_bound[i] <= snp_total) {         # If not yet the end of the VCF 
      
      fenetre = which(vec_snps_pos >= low_bound[i] & vec_snps_pos <= upper_bound[i])
      
      # ---- #### Convert the VCF to genind ...
      data_moment <- VCF1[fenetre,]
      genotype <- extract.gt(data_moment, element = "GT", mask = FALSE, as.numeric=F,
                             return.alleles = FALSE, IDtoRowNames = TRUE, 
                             extract = TRUE, convertNA = FALSE)
      
      geno_temp <- t(genotype)                      #### On transpose car objet genind a besoin d'une matrice 
      #### avec individus en lignes et SNPs en colonnes
      for (j in 1:ncol(geno_temp)){                 #### ./. --> NA
        geno_temp[geno_temp[,j]=="./.",j] <- NA 
      }
      mat_data <- df2genind(geno_temp, ploidy = 2, sep="/", pop = vectNom) 
      
      # ---- #### Compute the PCA ... 
      whale_pca <- indpca(mat_data)
      new_eigen <- whale_pca$ipca$eig / sum(whale_pca$ipca$eig)
      
      # SCAN VARIANCE : -----------
      scan_acp_var <- data.frame(PC1 = new_eigen[1],          # Récup % de variance des premiers axes de PC 
                                 PC2 = new_eigen[2],      
                                 low_f = low_bound[i],               
                                 upper_f = upper_bound[i], 
                                 Individu = vectNom, 
                                 Pop = names_ind$Population)
      final_acp_var1_2 <- rbind(final_acp_var1_2, scan_acp_var) 
      
      
      # SCAN POSITION : ------------
      scan_acp_pos <- data.frame(PC1 = whale_pca$ipca$li$Axis1,      # Récupérer les positions des individus
                             PC2 = whale_pca$ipca$li$Axis2,                     # sur les 2 premières PC
                             low_f = low_bound[i],               # Position du low_bound
                             upper_f = upper_bound[i],           # Position du upper_bound 
                             Individu = vectNom, 
                             Pop = names_ind$Population)
      
      final_acp_coord1_2 <- rbind(final_acp_coord1_2, scan_acp_pos) 
      
      
    }
  }
  # ---- #### Finalize the data.frames : -------------------------------------------------------
  scan_acp_var1_2 <- final_acp_var1_2 |>
    mutate(position = (low_f + upper_f)/2) |>
    pivot_longer(cols = c(PC1, PC2), 
                 names_to = "PC", 
                 values_to = "Variance")
  
  scan_acp_coord1_2 <- final_acp_coord1_2 |>
    mutate(position = (low_f + upper_f)/2) |>
    pivot_longer(cols = c(PC1, PC2), 
                 names_to = "PC", 
                 values_to = "coords")
  
  # ---- #### Save the data.frames : -----------------------------------------------------------
  
  saveRDS(scan_acp_var1_2, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/PCA/", chrom, "_scan_pca_var_pos.RDS", sep = ""))
  saveRDS(scan_acp_coord1_2, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/PCA/", chrom, "_scan_pca_coords.RDS", sep = ""))
  
  # ---- #### Plot it : ------------------------------------------------------------------------
  pdf(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/PCA/", chrom, "_scan_pca_var_pos.pdf", sep = ""))
  scan_acp_var1_2 |>
    ggplot(aes(x = position, y = Variance, color = PC)) + 
    geom_point(alpha = 0.3) +
    scale_y_continuous(limits = c(0, 0.5)) + 
    scale_color_manual(values = c("PC1" = "black", "PC2" = "red"))
  dev.off()
  
  pdf(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/PCA/", chrom, "_scan_pca_coord.pdf", sep = ""))
  scan_acp_coord1_2 |>
    filter(PC == "PC1") |>
    ggplot(aes(x = position, y = coords, color = Pop)) + 
    geom_line(alpha = 0.5) +
    #scale_color_manual(values = c("N_PACIFIC" = "darkgreen", "SOUTH_HE" = "red", "N_ATL" = "royalblue")) +
    theme_bw() + 
    labs(x = "Positions", 
         y = "PC1 coordinates", 
         title = paste("Scan PCA : 86 individuals PC1 coords along the chr", chrom))
  dev.off()
  
}
#########################################################################################################





