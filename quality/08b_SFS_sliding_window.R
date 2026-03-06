#################### 27/01/26 #####################
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

# Data : 
VCF1 <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered.RDS")
list_pop <- readRDS("data/list_pop.RDS")

source("quality/functions_for_td.r")

# ATTENTION : IMPORTER LES FONCTIONS NÉCESSAIRES 

################################ SLIDING WINDOW ################################

maf <- 0

#### Préparation des bornes à chaque itération : -------------------------------

slide <- 25000
window <- 100000

vec_pos <- getPOS(VCF1)
n_pop <- length(list_pop)
low_bound <- seq(1, max(vec_pos) + slide, slide)
upper_bound <- seq(window, max(vec_pos) + slide + window, slide)


#### Pour chaque fenêtre : -----------------------------------------------------

sfs_spectre_sliding <- c()        

for (i in 1:length(low_bound)) {
  fenetre = which(vec_pos > low_bound[i] & vec_pos < upper_bound[i])
  
  # ------- #### ... s'il y a des SNPs dans la fenêtre de positions ...  ----------
  if (length(fenetre) > 0) {
    data_temp <- VCF1[fenetre,]         # On sélectionne cette fenêtre de SNPs
    
    # ------- #### Et le SFS folded ... ---------    
    temp_spectre <- c()
    for (z in 1:length(list_pop)){
      data_bin <- vcfR2DNAbin(data_temp[, c("FORMAT",list_pop[[z]])], extract.indels = TRUE, 
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
positions_plot <- (upper_bound+low_bound)/2
sfs_spectre_sliding
sfs_spectre_final <- as.data.frame(na.omit(cbind(positions_plot, sfs_spectre_sliding)))
write_csv(sfs_spectre_final, "/shared/projects/multiwhaling/multiwhaling/plot/SFS/sfs_spectre_sliding_1e5_25000.csv")

# Calculer SFS norm par pop : 
a <- 0
list_plot_div <- list()
div_gen_final <- c()

for (i in seq_along(list_pop)) {
  
  # ----------  #### Sélection de la pop : 
  n_cols <- length(list_pop[[i]])      # Nb d'individus appartenant à la pop i
  col_range <- (a + 1):(a + n_cols)    # Sélectionne les individus (colonnes) à partir de la colonne suivant la dernière de la pop précédente ... dernier de la pop
  
  # ----------  #### Calcul des indices de diversité : TD, Pi, W et S
  matrice_sfs_bypop <- sfs_spectre_sliding[, col_range]
  qq = apply(matrice_sfs_bypop, 1,calcola_TD_folded)          # Calcul du SFS ligne par ligne          
  divgen = matrix(unlist(qq), ncol= 4, nrow = nrow(matrice_sfs_bypop), byrow = T)
  divgen_spectre <- as_tibble(cbind(positions_plot, divgen))
  colnames(divgen_spectre) <- c("positions", "TD", "Pi", "W", "S")
  
  
  # ---------- #### Faire le scaling par le nombre de positions comprises dans chaque fenêtre 
  divgen_spectre <- divgen_spectre |>
    mutate(Pi = Pi/window,            # Pour les deux estimateurs, on scale par le nb total de positions séquencées (SNP & NPP)
           W = W/window,              # Ici, ça correspond à la taille de la fenêtre sélectionnée
           Pop = paste(names(list_pop[i])))              
  
  saveRDS(divgen_spectre, paste("/shared/projects/multiwhaling/multiwhaling/plot/SFS/sliding_SFS_1e5_25000_", names(list_pop[i]), ".RDS", sep = ""))
  div_gen_final <- rbind(div_gen_final, divgen_spectre)
  
  # ----------  #### Plot des indices le long du chromosome : 
  plot_div <- divgen_spectre |>
    pivot_longer(cols = 2:5, 
                 names_to = "indices_div", 
                 values_to = "value") |>
    ggplot(aes(x = positions, y = value, color = indices_div)) +
    geom_point() +
    facet_wrap(~ indices_div, scales = "free_y") +
    labs(title = paste(names(list_pop[i])), 
         x = "Positions", 
         y = "Value", 
         color = "Diversity indices")
  list_plot_div[[i]] <- plot_div
  
  # ---------- #### On passe à la pop suivante : 
  a <- a + n_cols
}

saveRDS(div_gen_final, "/shared/projects/multiwhaling/multiwhaling/plot/SFS/div_gen_final_1e5_25000.RDS")
saveRDS(list_plot_div, "/shared/projects/multiwhaling/multiwhaling/plot/SFS/list_plot_div_1e5_25000.RDS")
