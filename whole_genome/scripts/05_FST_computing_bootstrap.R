############################# 20/03/26 ##################################
############### Compute pairwise FST and permutations ###################
######################## Achille RIVOLTELLA #############################

# Library : 
library(vcfR)
library(rlist)
library(ggplot2)
library(tidyverse)
library(pheatmap)

#######################################################################################################
# ARGUMENTS à mettre en entrée :                                                                      #
#         - VCF avec SNPs ;                                                                           #
#         - objet (liste) contenant les noms des pops et les individus associés                       #
#         - n° chromosome                                                                             #
#         - Nombre de permutations voulues                                                            #  
#                                                                                                     #
# SORTIES de la fonction :                                                                            #
#         - Fst par paires de pop observés et plot associé                                            #
#         - p_values pour chaque paire, obtenue par permutation des pop                               #
#         - Table des fst pour chaque permutation                                                     #
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

permut = args[2]

### Fonctions sources :
source("/shared/projects/multiwhaling/multiwhaling/whole_genome/fonctions/fonctions.R")

fst_pairwise_bootstrap <- function(VCF1, list_pop, names_ind, permut) {  # -----------------------------# 
  
  #------------------------------------------------------------------------------------------------------
  ############################### FST observé des populations deux à deux ###############################
  #------------------------------------------------------------------------------------------------------
  
  # Préparation des objets nécessaires : 
  pairs <- combn(seq_along(list_pop), 2)  # Création des paires de pop pour lesquelles on calcule le FST
  
  fst_moy <- data.frame(matrix(NA, nrow = length(list_pop), ncol = length(list_pop), 
                               dimnames = list(names(list_pop), names(list_pop))))
  
  tab_fst_obs <- data.frame(matrix(NA, nrow = ncol(pairs), ncol = 2))
  colnames(tab_fst_obs) <- c("Combi", "FST")
  
  # Calcul du FST observé par paires : 
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]
    j <- pairs[2, k]
    
    list_pop_hudson <- list(list_pop[[i]], list_pop[[j]])
    fst_pairwise <- fst.hudson(VCF1, list_pop_hudson)
    
    fst_moy[i, j] <- unlist(fst_pairwise[1])         # Valeur moyenne de fst entre 2 pop
    
    tab_fst_obs[k, 1] <- paste(i, j, sep = "_")
    tab_fst_obs[k, 2] <- unlist(fst_pairwise[1]) 
  }
  write.csv(fst_moy, file = "/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/05_FST/fst_obs.csv")
  write.csv(tab_fst_obs, file = "/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/05_FST/fst_obs_long.csv")
  
  png("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/05_FST/fst_pairwise_heatmap.png", height = 750, width = 750)
  pheatmap(fst_moy, cluster_rows=F, cluster_cols=F, na_col="white",main = "Pairwise Fst",
           color = colorRampPalette(c("seashell1", "yellow", "firebrick3"))(50))
  dev.off()
  
  
  
  #-----------------------------------------------------------------------------------------------------#
  ###################################### Permutations sur les FST #######################################
  #-----------------------------------------------------------------------------------------------------#
  
  # Nombre de permutations à effectuer : 
  print(permut)
  
  # On crée le dataframe qui va accueillir les données de FST à chaque itération de bootstrap : 
  tab_fst_comb <- data.frame(matrix(NA, nrow = ncol(pairs), ncol = 2))
  colnames(tab_fst_comb) <- c("Combi", "FST")
  
  # Dataframe final issu de la fusion des valeurs de FST à chaque itération de bootstrap : 
  tab_fst_all <- c()
  
  for (n in 1:permut) {
    # Permutter les individus entre les pops : 
    names_ind_perm <- names_ind     # Pour garder la liste originale pour faire le sample à chaque permutation 
    names_ind_perm$Individu <- sample(names_ind$Individu)
    
    # Liste de pop mélangées : 
    list_bootstrap <- split(names_ind_perm$Individu, names_ind_perm$Population)    
    
    # Calcul du FST par paires : 
    for (k in seq_len(ncol(pairs))) {
      i <- pairs[1, k]
      j <- pairs[2, k]
      
      list_pop_hudson <- list(list_bootstrap[[i]], list_bootstrap[[j]])
      fst_pairwise <- fst.hudson(VCF1, list_pop_hudson)
      
      tab_fst_comb[k, 1] <- paste(i, j, sep = "_")
      tab_fst_comb[k, 2] <- unlist(fst_pairwise[1])       # 1 ligne = 1 combinaison de 2 pops et 1 FST :
    }
    tab_fst_all <- rbind(tab_fst_all, tab_fst_comb)       
  }
  write.csv(tab_fst_all, file = "/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/05_FST/fst_permut.csv")
  
  
  #-----------------------------------------------------------------------------------------------------#
  ################################### Calcul des IC95 et p_value : ######################################
  #-----------------------------------------------------------------------------------------------------#
  
  tab_fst_obs                                 # FST observés 
  tab_stat <- tab_fst_all |>                  # FST permutations (si pop sont panmictiques)
    group_by(Combi) |>
    reframe(mean_cl_boot(FST), 
            n_perm = n(),
            b = sum(FST >= tab_fst_obs$FST[match(Combi, tab_fst_obs$Combi)]),  # b = nombre de permutations avec FST_permut ≥ FST_obs
            p_value = (b + 1) / (n_perm + 1))
  
  #### Création du tableau des p_value (même forme que fst_moy) : ----------------
  fst_pval <- data.frame(matrix(NA, nrow = length(list_pop), ncol = length(list_pop), 
                                dimnames = list(names(list_pop), names(list_pop))))
  n_combi <- seq_along(tab_fst_comb$Combi)
  
  for (p in n_combi) {
    # On extrait l'index des 2 pops 
    i <- as.numeric(str_split(tab_fst_comb$Combi[p], "_")[[1]][1])
    j <- as.numeric(str_split(tab_fst_comb$Combi[p], "_")[[1]][2])
    
    # Et on assigne la p_value correspondant dans un df de la même forme que fst_moy : 
    fst_pval[i, j] <- tab_stat$p_value[p]
  }
  write.csv(fst_pval, file = "/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/05_FST/fst_pval.csv")
}
###############################################################################################################

fst_pairwise_bootstrap(VCF, list_pop, names_ind, permut)



