################### 20/03/26 ######################
################# PCA analysis ####################
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)
library(pegas)
library(LEA)

#######################################################################################################
############################## Version avec seulement qq pops #########################################
#######################################################################################################

#######################################################################################################
# ARGUMENTS à mettre en entrée :                                                                      #
#         - VCF avec SNPs ; VCF avec toutes les positions ;                                           #
#         - objets (liste et tibble) contenant les noms des pops et les individus associés            #
#         - n° chromosome                                                                             #
#                                                                                                     #
# SORTIES de la fonction :                                                                            #
#         - VCF transformé en .geno                                                                   #
#         - Dossier ACP avec eigenvalues, projections sur les axes                                    #
#         - Plot des eigenvalues et plot PC1 & PC2                                                    #
#         - Datasets d'ancestry par valeurs de K                                                      #
#         - Plots de structure                                                                        #
#         - Projets de PCA et sNMF                                                                    # 
#######################################################################################################


#### ARGUMENTS :  
args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]
#chrom <- 21

VCF = read.vcfR(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/VCF_filtered_maf_", 
                      chrom, "_LD_Pruned.vcf.gz", sep = ""))   
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/list_pop.RDS")

names_ind <- readRDS("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/names_ind.RDS")


pca_clustering <- function(VCF1, list_pop, names_ind){
  
  print(ocean)
  
  #------------------------------------------------------------------------------------------------------
  #################################### PRINCIPAL COMPONENT ANALYSIS #####################################
  #------------------------------------------------------------------------------------------------------
  
  #### Transformer le fichier VCF en .GENO : ------------------------------------------------------------
  fun_geno_allele<-function(data){
    for (i in 1:length(data)) {
      if (data[i]=="./.") {data[i]<-"9"}
      if (data[i]=="0/0") {data[i]<-0}
      if (data[i]=="0/1") {data[i]<-1}
      if (data[i]=="1/1") {data[i]<-2}
    }
    return(data)
  }
  
  genotype <- extract.gt(VCF1, element = "GT", mask = FALSE, as.numeric=F, return.alleles = FALSE, 
                         IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) 
  n_ind <- ncol(genotype)
  snps_tot <- nrow(genotype)
  qqq <- apply(genotype, 1, fun_geno_allele)
  genotype_num <- matrix(as.numeric(qqq), ncol = n_ind, nrow = snps_tot, byrow = T)
  
  write.table(genotype_num, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/VCF_", chrom, "_", ocean, ".geno", sep = ""), 
              quote = F, sep = "", row.names = F, col.names=F)
  
  #### PCA : --------------------------------------------------------------------------------------------
  whale_pca <- pca(input.file = paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/VCF_", chrom, "_", ocean, ".geno", sep = ""))
  #summary(whale_pca)
  print("PCA : Running PCA")
  
  projectpca <- load.pcaProject(paste("VCF_", chrom,"_", ocean, ".pcaProject", sep = ""))
  
  # Eigenvalues : 
  png(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/eigenvalues_", chrom, "_", ocean,".png", sep = ""))
  plot(whale_pca, lwd = 5, col = "blue", cex = .7, xlab = ("Factors"), ylab = "Eigenvalues")
  dev.off()
  print("PCA : Printed eigenvalues")
  
  #### MISE EN FORME DES RÉSULTATS : -------------------------------------------------------------------
  
  # On récupère les données dans un nouveau dossier généré automatiquement par la fonction 'pca' utilisée
  eigen_pca <- read.table(paste(projectpca@projDir, projectpca@pcaDir, projectpca@eigenvalue.file, sep = ""))
  eigen_pca <- (eigen_pca / ((sum(eigen_pca)))) * 100
  eigen_pca <- unlist(round(eigen_pca, digits = 2))
  saveRDS(eigen_pca, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/eigenvalues_", chrom, "_", ocean, ".RDS", sep = ""))
  
  #### SAUVEGARDE DES PLOTS : --------------------------------------------------------------------------
  
  # Récupération des coordonnées de chaque individu dans les projections : 
  proj <- read.table(paste(projectpca@projDir, projectpca@pcaDir, projectpca@projection.file, sep = ""), header = F)
  proj <- as_tibble(proj) |>
    dplyr::select(1:4) |>
    mutate(pop = names_ind |> arrange(Population))
  saveRDS(proj, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/projections_ACP_", chrom, "_", ocean, ".RDS", sep = ""))
  
  png(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/ACP_", chrom, "_", ocean, ".png", sep = ""))
  proj |>
    ggplot(aes(x = V1, y = V2, color = pop$Population)) + 
    geom_point(size = 2) + 
    theme_bw() + 
    labs(x = paste("PC1 ", as.character(eigen_pca[1]), "%"), 
         y = paste("PC2 ", as.character(eigen_pca[2]), "%"), 
         color = "Population")
  dev.off()
  print("PCA : Saved plots")
  
  #------------------------------------------------------------------------------------------------------
  ############################################# sNMF CLUSTERING #########################################
  #------------------------------------------------------------------------------------------------------
  
  #### CLUSTERING : -------------------------------------------------------------------------------------
  whale_snmf <- snmf(input.file= paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/VCF_", chrom, "_", ocean, ".geno", sep = ""), 
                     K= 1:length(list_pop), repetitions = 3, project="new", entropy=T)
  print("sNMF : Running sNMF")
  
  projectsNMF <- load.snmfProject(file = paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/VCF_", chrom, "_", ocean, ".snmfProject", sep = ""))
  
  png(file = paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/entropy_K_", chrom, "_", ocean, ".png", sep = ""))
  plot(whale_snmf, col = "blue", pch = 19, cex = 1.2) ### Identifier le K avec l'entropie la plus faible
  dev.off()
  print("sNMF : Printed entropy")
  
  #### MISE EN FORME DES RÉSULTATS : -------------------------------------------------------------------
  list_plot <- c()
  adm_coeff_all <- c()
  order <- c("BER", "KAR", "CHI", "PER", "MAD", "WAP", "SPM", "GOM")
  
  print("Plotting results by number of K")
  for (i in 2:length(list_pop)) {
    best <- which.min(cross.entropy(whale_snmf, K = i))
    adm_coeff <- Q(whale_snmf, K = i, run = best)
    
    adm_coeff <- as_tibble(adm_coeff) |>
      mutate(ind = names_ind$Individu, 
             pop = factor(names_ind$Population, levels = order)) |>
      arrange(pop) |>
      pivot_longer(cols = 1:i, names_to = "K", values_to = "Ancestry_prop") |>
      mutate(K = gsub("V", "", K), 
             K_i = i)
    
    adm_coeff_all <- rbind(adm_coeff_all, adm_coeff)
    
    plot_structure <- adm_coeff |>
      mutate(ind = factor(ind, levels = unique(adm_coeff$ind))) |>
      ggplot(aes(x = ind, y = Ancestry_prop, fill = K, color = K)) + 
      geom_col() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3, color = "black")) + 
      labs(x = "Individu", 
           y = "Ancestry proportions")
    
    list_plot <- c(list_plot, list(plot_structure))
  }
  
  saveRDS(list_plot, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/by_pop/plot_struct_", chrom, "_", ocean, ".RDS", sep = ""))
  saveRDS(adm_coeff_all, paste(projectsNMF@projDir, projectsNMF@snmfDir, "admix_ancestry_", chrom, "_", ocean, ".RDS", sep = ""))
  print("sNMF : Saved plots and data")
}
#####################################################################################################

# # Pour Nord Atlantique : 
# list_pop <- list_pop[6:7]
# names_ind <- enframe(list_pop, name = "Population", value = "Individu") |>
#   unnest(Individu)
# 
# VCF <- VCF[,c("FORMAT", unlist(list_pop))] 
# # Pour hémisphère Sud :
# list_pop <- list_pop[c(3:5, 8)]
# names_ind <- enframe(list_pop, name = "Population", value = "Individu") |>
#   unnest(Individu)
# VCF <- VCF[, c("FORMAT", unlist(list_pop))]

# Assigner les bassins océaniques à chaque population : 
names_ind_basins <- names_ind |>
  mutate(basin = ifelse(Population == "BER" | Population == "KAR", "NORTH_PACIFIC", 
                      ifelse(Population == "PER" | Population == "CHI" | Population == "MAD" | Population == "WPA", "SOUTH_HEMISPHERE", "NORTH_ATLANTIC")))
list_pop_basins <- split(names_ind_basins$Individu, names_ind_basins$basin)
basins <- names(list_pop_basins)

for (ocean in basins) {
  # Filtrer par ocean : 
  names_ind_tmp <- names_ind_basins |>
    filter(basin == ocean) |>
    mutate(Population = as.character(Population))
  list_pop_tmp <- split(names_ind_tmp$Individu, names_ind_tmp$Population)

  # Appliquer le filtre sur le VCF :
  colnames(VCF@gt)
  VCF <- VCF[, c("FORMAT", unlist(list_pop_tmp))]
  colnames(VCF@gt)

  # Exécuter le code à chaque fois :
  pca_clustering(VCF, list_pop, names_ind_tmp)
}


