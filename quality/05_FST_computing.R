############################# 27/01/26 ##################################
############### Compute pairwise FST and permutations ###################
######################## Achille RIVOLTELLA #############################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)
library(pegas)
library(pheatmap)


#### Read the data :  
VCF1 <- read.vcfR("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_maf_LD_Pruned.vcf.gz")
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/list_pop.RDS")
names_ind <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/names_ind.RDS")

source("/shared/projects/multiwhaling/multiwhaling/quality/functions_for_td.r")

#------------------------------------------------------------------------------------------------------
########## FST observé des populations deux à deux ########## -----------------------------------------
#------------------------------------------------------------------------------------------------------

# Data nécessaires : 
pairs <- combn(seq_along(list_pop), 2)  # Création des paires de pop pour lesquelles on calcule le FST

fst_moy <- data.frame(matrix(NA, nrow = length(list_pop), ncol = length(list_pop), 
                             dimnames = list(names(list_pop), names(list_pop))))

tab_fst_obs <- data.frame(matrix(NA, nrow = ncol(pairs), ncol = 2))
colnames(tab_fst_obs) <- c("Combi", "FST")


for (k in seq_len(ncol(pairs))) {
  i <- pairs[1, k]
  j <- pairs[2, k]
  
  list_pop_hudson <- list(list_pop[[i]], list_pop[[j]])
  fst_pairwise <- fst.hudson(VCF1, list_pop_hudson)
  
  fst_moy[i, j] <- unlist(fst_pairwise[1])                      # Valeur moyenne de fst entre 2 pop
  
  tab_fst_obs[k, 1] <- paste(i, j, sep = "_")
  tab_fst_obs[k, 2] <- unlist(fst_pairwise[1]) 
 }

write.csv(fst_moy, file = "/shared/projects/multiwhaling/multiwhaling/plot/FST/mean_fst_pop.csv")

png("/shared/projects/multiwhaling/multiwhaling/plot/FST/fst_pairwise_heatmap.png", height = 750, width = 750)
pheatmap(fst_moy, cluster_rows=F, cluster_cols=F, na_col="white",main = "Pairwise Fst",
         color = colorRampPalette(c("seashell1", "yellow", "firebrick3"))(50))
dev.off()



#------------------------------------------------------------------------------------------------------
########## Permutations sur les FST ########## --------------------------------------------------------
#------------------------------------------------------------------------------------------------------

# Nombre de permutations à effectuer : 
permut <- 1:1000

# On crée le dataframe qui va accueillir les données de FST à chaque itération de bootstrap : 
tab_fst_comb <- data.frame(matrix(NA, nrow = ncol(pairs), ncol = 2))
colnames(tab_fst_comb) <- c("Combi", "FST")

# Dataframe final issu de la fusion des valeurs de FST à chaque itération de bootstrap : 
tab_fst_all <- c()

for (n in permut) {
  # Permutter les individus entre les pops : 
  names_ind_perm <- names_ind     # Pour garder la liste originale pour faire le sample à chaque permutation 
  names_ind_perm$value <- sample(names_ind$value)
  
  # Liste de pop mélangées : 
  list_bootstrap <- split(names_ind_perm$value, names_ind_perm$pop)    
  
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

#### Calcul des IC95 et p_value : ----------------------------------------------

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


#### Sauvegarder des données obtenues : ----------------------------------------
saveRDS(tab_fst_obs, "/shared/projects/multiwhaling/multiwhaling/plot/FST/fst_obs.RDS")
saveRDS(tab_fst_all, "/shared/projects/multiwhaling/multiwhaling/plot/FST/fst_permut.RDS")
saveRDS(tab_stat, "/shared/projects/multiwhaling/multiwhaling/plot/FST/fst_stats.RDS")
write.csv(fst_pval, "/shared/projects/multiwhaling/multiwhaling/plot/FST/fst_pval.csv")

