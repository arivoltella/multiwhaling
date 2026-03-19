#################### 11/03/26 #####################
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
VectPosFiltered <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/vect_position/pos_all.RDS")         # Toutes positions du chromosomes 
list_pop <- readRDS("data/list_pop.RDS")

source("quality/functions_for_td.r")

# ATTENTION : IMPORTER LES FONCTIONS NÉCESSAIRES 
sliding_spectre_onepop<-function(dati, window, slide){
  dim_locus<-window
  interval<-slide
  vettore_snps_pos<-getPOS(dati)
  vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
  vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
  n_ind<-length(vettore_nomi)
  low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
  upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
  per_plot<-(upper_bound+low_bound)/2
  
  ris_spectre_sliding<-c()
  for (i in 1:length(low_bound)) {
    sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
    if (length(sss)>0) {
      dati_temp<-dati[sss,]
      temp_spectre<-c()
      dati_bin<-vcfR2DNAbin(dati_temp, extract.indels = TRUE, consensus = FALSE,  unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL, verbose = TRUE) ### trasformo in pratica il vcf in sequenze, e posso quindi calcolare sia il SFS che la nucleotide diversity. Ma non l'heterozygosity
      sfs_tot<-site.spectrum(dati_bin, folded = TRUE) #### sfs folded con funzione PEGAS
      temp_spectre<-c(temp_spectre,sfs_tot)
      ris_spectre_sliding<-rbind(ris_spectre_sliding,temp_spectre)
    } else {
      ris_spectre_sliding<-rbind(ris_spectre_sliding, rep(NaN,n_ind))
    }
  }
  
  finale_spectre<-cbind(per_plot,ris_spectre_sliding)
  finale_spectre<-na.omit(finale_spectre)
  
  return(finale_spectre)
}

#### FILTRER LES NA POUR ANALYSE SFS ####

################################ SLIDING WINDOW ################################
z = 1 # Bering

# Filtre pour 1 pop : 
data <- VCF1[, c("FORMAT",list_pop[[z]])]

# sliding window sur 1 pop 
test_100K <- sliding_spectre_onepop(data, 500000, 10000)
test_100K <- as.data.frame(test_100K)
write_csv(test_100K, "/shared/projects/multiwhaling/multiwhaling/plot/scan_genom/SFS/test_100K.csv")

# test_50K <- sliding_spectre_onepop(data, 50000, 10000)
# test_50K <- as.data.frame(test_50K)
# write_csv(test_50K, "/shared/projects/multiwhaling/multiwhaling/plot/scan_genom/SFS/test_50K.csv")
# 
# 
# 
# ######################################################
# ######### Suite script SFS-sliding_fct_S #############
# ######################################################
# 
# # ----------  #### Calcul des indices de diversité : TD, Pi, W et S
# test_100K <- read_csv("~/Desktop/scan_genom/SFS/test_100K.csv")
# test_50K <- read_csv("~/Desktop/scan_genom/SFS/test_50K.csv")
# 
# positions_plot_100K <- test_100K[,1]
# positions_plot_50K <- test_50K[,1]
# 
# test_100K <- test_100K[,-1]
# test_50K <- test_50K[,-1]
# 
# # Récupération du nb de SNP moyen par fenêtre : 
# qq = apply(test_100K, 1,calcola_TD_folded)                  
# divgen = matrix(unlist(qq), ncol= 4, nrow = nrow(test_100K), byrow = T)
# divgen_spectre_100K <- as_tibble(cbind(positions_plot_100K, divgen))
# colnames(divgen_spectre_100K) <- c("positions", "TD", "Pi", "W", "S")
# 
# mean(divgen_spectre_100K$S)
# 
# qq = apply(test_50K, 1,calcola_TD_folded)          # Calcul du SFS ligne par ligne          
# divgen = matrix(unlist(qq), ncol= 4, nrow = nrow(test_50K), byrow = T)
# divgen_spectre_50K <- as_tibble(cbind(positions_plot_50K, divgen))
# colnames(divgen_spectre_50K) <- c("positions", "TD", "Pi", "W", "S")
# 
# mean(divgen_spectre_50K$S)
# 
# 
# # ------------------------------------------------------------------------------
# ########## Calcul distance SFS fenêtre et SFS total, par population ############
# # ------------------------------------------------------------------------------
# 
# # Pour Bering : 
# norm_sfs_bypop <- c(0.2601204, 0.2778992, 0.2899138, 0.2662977, 0.2760875, 0.2928641,
#                     0.2884587, 0.3032762, 0.2689080, 0.2542088, 0.2819353)
# 
# ################################### 100K ########################################
# 
# n_cols <- length(list_pop[[z]])      # Nb d'individus appartenant à la pop i
# 
# # ----------  #### Calcul de la distance euclidienne entre SFS norm total et par fenêtre : 
# eeeeee_norep <- t(apply(test_100K, 1, calcola_normalized_foldedSFS))
# sfsnorm_bywindow <- eeeeee_norep[, -c(1:n_cols)]
# 
# dist_euclid <- apply(sfsnorm_bywindow, 1, calcola_dist_spettro, unlist(norm_sfs_bypop))
# 
# dist_euclid_bypop_100K <- c(list(dist_euclid))
# 
# # Récupération du nb de SNP par fenêtre : 
# qq = apply(test_100K, 1,calcola_TD_folded)          # Calcul du SFS ligne par ligne          
# divgen = matrix(unlist(qq), ncol= 4, nrow = nrow(test_100K), byrow = T)
# divgen_spectre <- as_tibble(cbind(positions_plot_100K, divgen))
# colnames(divgen_spectre) <- c("positions", "TD", "Pi", "W", "S")
# n_snp <- divgen_spectre$S
# mean((n_snp))
# 
# 
# dist_euclid_bypop_100K[[2]] <- positions_plot_100K
# final_dist_euclid_100K <- bind_cols(dist_euclid_bypop_100K)
# final_dist_euclid_100K |>
#   ggplot(aes(x = per_plot, y = ...1)) + 
#   geom_line()
# 
# ################################### 50K ########################################
# 
# n_cols <- length(list_pop[[z]])      # Nb d'individus appartenant à la pop i
# 
# # ----------  #### Calcul de la distance euclidienne entre SFS norm total et par fenêtre : 
# eeeeee_norep <- t(apply(test_50K, 1, calcola_normalized_foldedSFS))
# sfsnorm_bywindow <- eeeeee_norep[, -c(1:n_cols)]
# 
# dist_euclid <- apply(sfsnorm_bywindow, 1, calcola_dist_spettro, unlist(norm_sfs_bypop))
# 
# dist_euclid_bypop_50K <- c(list(dist_euclid))
# 
# # Récupération du nb de SNP par fenêtre : 
# qq = apply(test_50K, 1,calcola_TD_folded)          # Calcul du SFS ligne par ligne          
# divgen = matrix(unlist(qq), ncol= 4, nrow = nrow(test_50K), byrow = T)
# divgen_spectre <- as_tibble(cbind(positions_plot_50K, divgen))
# colnames(divgen_spectre) <- c("positions", "TD", "Pi", "W", "S")
# n_snp <- divgen_spectre$S
# mean(n_snp)
# 
# dist_euclid_bypop_50K[[2]] <- positions_plot_50K
# final_test_dist_50K <- bind_cols(dist_euclid_bypop_50K)
# final_test_dist_50K |>
#   ggplot(aes(x = per_plot, y = ...1)) + 
#   geom_line()
