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
#chrom <- args[3]
chrom <- 21

VCF = read.vcfR(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/VCF_filtered_", 
                      chrom, ".vcf.gz", sep = ""))   
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/list_pop.RDS")
names_ind <- readRDS("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/names_ind.RDS")


pca_clustering <- function(VCF1, list_pop, names_ind){
  
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
  
  write.table(genotype_num, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/geno/VCF_", chrom, ".geno", sep = ""), 
              quote = F, sep = "", row.names = F, col.names=F)
  
  #### PCA : --------------------------------------------------------------------------------------------
  whale_pca <- pca(input.file = paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/geno/VCF_", chrom, ".geno", sep = ""))
  summary(whale_pca)
  
  projectpca <- load.pcaProject(paste("VCF_", chrom, ".pcaProject", sep = ""))
  
  # Eigenvalues : 
  pdf(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/eigenvalues_", chrom, ".pdf", sep = ""))
  plot(whale_pca, lwd = 5, col = "blue", cex = .7, xlab = ("Factors"), ylab = "Eigenvalues")
  dev.off()
  
  
  #### MISE EN FORME DES RÉSULTATS : -------------------------------------------------------------------
  
  # On récupère les données dans un nouveau dossier généré automatiquement par la fonction 'pca' utilisée
  eigen_pca <- read.table(paste(projectpca@projDir, projectpca@pcaDir, projectpca@eigenvalue.file, sep = ""))
  eigen_pca <- (eigen_pca / ((sum(eigen_pca)))) * 100
  eigen_pca <- unlist(round(eigen_pca, digits = 2))
  
  
  #### SAUVEGARDE DES PLOTS : --------------------------------------------------------------------------
  
  # Récupération des coordonnées de chaque individu dans les projections : 
  proj <- read.table(paste(projectpca@projDir, projectpca@pcaDir, projectpca@projection.file, sep = ""), header = F)
  proj <- as_tibble(proj) |>
    dplyr::select(1:4) |>
    mutate(pop = names_ind |> arrange(Population))
  saveRDS(proj, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/projections_ACP_", chrom, ".RDS", sep = ""))
  
  pdf(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/ACP_", chrom, ".pdf", sep = ""))
  proj |>
    ggplot(aes(x = V1, y = V2, color = pop$Population)) + 
    geom_point(size = 2) + 
    theme_bw() + 
    labs(x = paste("PC1 ", as.character(eigen_pca[1]), "%"), 
         y = paste("PC2 ", as.character(eigen_pca[2]), "%"), 
         color = "Population")
  dev.off()
  
  
  #------------------------------------------------------------------------------------------------------
  ############################################# sNMF CLUSTERING #########################################
  #------------------------------------------------------------------------------------------------------
  
  #### CLUSTERING : -------------------------------------------------------------------------------------
  whale_snmf <- snmf(input.file= "/shared/projects/multiwhaling/multiwhaling/whole_genome/data/geno/VCF_21.geno", 
                     K= 1:6, repetitions = 3, project="new", entropy=T)
  summary(whale_snmf)
  projectsNMF <- load.snmfProject(file = paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/geno/VCF_", chrom, ".snmfProject", sep = ""))
  
  
  pdf(file = "/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/entropy_K_", chrom, ".pdf")
  plot(whale_snmf, col = "blue", pch = 19, cex = 1.2) ### Identifier le K avec l'entropie la plus faible
  dev.off()
  
  #### MISE EN FORME DES RÉSULTATS : -------------------------------------------------------------------
  list_plot <- c()
  order <- c("BERING", "KARAGINSKY", "CHILI", "PEROU", "MADAGASCAR", "N_ATL")
  
  for (i in 2:6) {
    best <- which.min(cross.entropy(whale_snmf, K = i))
    adm_coeff <- Q(whale_snmf, K = i, run = best)
    
    adm_coeff <- as_tibble(adm_coeff) |>
      mutate(ind = names_ind$Individu, 
             pop = factor(names_ind$Population, levels = order)) |>
      arrange(pop) |>
      pivot_longer(cols = 1:i, names_to = "K", values_to = "Ancestry_prop") |>
      mutate(K = gsub("V", "", K))
    
    saveRDS(adm_coeff, paste(projectsNMF@projDir, projectsNMF@snmfDir, "ancestry_", chrom, "_K_", i, ".RDS", sep = ""))
    
    plot_structure <- adm_coeff |>
      mutate(ind = factor(ind, levels = unique(adm_coeff$ind))) |>
      ggplot(aes(x = ind, y = Ancestry_prop, fill = K, color = K)) + 
      geom_col() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3, color = "black")) + 
      labs(x = "Individu", 
           y = "Ancestry proportions")
    
    list_plot <- c(list_plot, list(plot_structure))
  }
  
  structure <- wrap_plots(list_plot) + plot_layout(axis_titles = "collect")
  saveRDS(structure, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/06_PCA_clustering/structure_", chrom, ".RDS"))
  
}
