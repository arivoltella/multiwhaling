#################### 16/03/26 #####################
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
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/list_pop.RDS")
names_ind <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/names_ind.RDS")

# ATTENTION : IMPORTER LES FONCTIONS NÉCESSAIRES 
source("/shared/projects/multiwhaling/multiwhaling/quality/functions_for_td.r")


# ------------------------------------------------------------------------------------------------------
############################### Calculs distance SFS et SFS par window #################################
# ------------------------------------------------------------------------------------------------------

# Mise en forme des données : 
sfs_spectre_sliding <- read_csv("/shared/projects/multiwhaling/multiwhaling/data/scan_genom/SFS/sfs_spectre_sliding_50_1e4.csv")
positions_plot <- unlist(sfs_spectre_sliding[,1])
sfs_spectre_sliding <- sfs_spectre_sliding[,-1]

# Récupérer le SFS total pour chaque pop : 
norm_sfs_bypop <- c()

# ------------------------------------------------------------------------------
############ Création des SFS normalized population par population #############
# ------------------------------------------------------------------------------
for (i in seq_along(list_pop)) {
  
  # ---- ####  Tri par population : 
  pop <- list_pop[i]
  data <- VCF1[,c("FORMAT", unlist(pop))] 

  # ---- #### Filtre des NAs : 
  genotypes <- extract.gt(data, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  NAs <- rowSums(genotypes == "./.")
  data <- subset(data, NAs < 1)
  
  # ---- #### Calcul du SFS : 
  dati_bin <- vcfR2DNAbin(data, extract.indels = TRUE, consensus = FALSE, 
                          unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL,
                          verbose = TRUE)
  sfs_tot <- as.vector(site.spectrum(dati_bin, folded = TRUE)) 
  
  # ---- #### Normalized SFS + modèle neutre :
  norm_sfs_tot <- calcola_normalized_foldedSFS(sfs_tot) 
  norm_sfs_bypop <- c(norm_sfs_bypop, list(norm_sfs_tot[,2]))
}
# Pour chaque pop, autant de valeurs de SFSnorm que d'individus 

# ------------------------------------------------------------------------------
########## Calcul distance SFS fenêtre et SFS total, par population ############
# ------------------------------------------------------------------------------

# Fonction de calcul distance euclidienne : 
calcola_dist_spettro <- function(data, vettore_average){
  distanza <- sum((data - vettore_average)^2)
  return(distanza)
}

# Calculer SFS norm par pop : 
a <- 0
dist_euclid_bypop <- c()

# ---- #### Pour chaque population : -------------------------------------------
for (z in seq_along(list_pop)) {
  
  # ----------  #### Sélection de la pop : 
  n_cols <- length(list_pop[[z]])      # Nb d'individus appartenant à la pop i
  col_range <- (a + 1):(a + n_cols)    # Sélectionne les individus (colonnes) à partir de la colonne suivant la dernière de la pop précédente ... dernier de la pop
  
  # ----------  #### Calcul de la distance euclidienne entre SFS norm total et par fenêtre : 
  matrice_sfs_bypop <- sfs_spectre_sliding[, col_range]
  eeeeee_norep <- t(apply(matrice_sfs_bypop, 1, calcola_normalized_foldedSFS))
  sfsnorm_bywindow <- eeeeee_norep[, -c(1:n_cols)]
  
  dist_euclid <- apply(sfsnorm_bywindow, 1, calcola_dist_spettro, unlist(norm_sfs_bypop[z]))
  
  dist_euclid_bypop <- c(dist_euclid_bypop, list(dist_euclid))
  
  # ---------- #### On passe à la pop suivante : -------------------------------
  a <- a + n_cols
}


#### Arrangement du tableau et plots : -----------------------------------------
dist_euclid_bypop[[7]] <- positions_plot
final_dist_euclid_bypop <- bind_cols(dist_euclid_bypop)
colnames(final_dist_euclid_bypop) <- c(unique(names_ind$Population), "position")

saveRDS(final_dist_euclid_bypop, "/shared/projects/multiwhaling/multiwhaling/plot/scan_genom/SFS/dist_euclid_SFS.RDS")

# final_dist_euclid_bypop |>
#   pivot_longer(cols = 1:6, 
#                names_to = "Population", 
#                values_to = "Eta_2") |> 
#   ggplot(aes(x = position, y = Eta_2, color = Population)) + 
#   geom_line() + 
# 
# 
# final_dist_euclid_bypop |>
#   pivot_longer(cols = 1:6, 
#                names_to = "Population", 
#                values_to = "Eta_2") |>
#   reframe(mean_eta = mean(Eta_2), 
#           .by = Population)

