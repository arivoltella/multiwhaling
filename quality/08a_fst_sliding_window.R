#################### 27/01/26 #####################
############ Sliding window FST : #################
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
VCF1 <- read.vcfR("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_maf_LD_Pruned.vcf.gz")
list_pop <- readRDS("data/list_pop.RDS")

source("quality/functions_for_td.r")

# ATTENTION : IMPORTER LES FONCTIONS NÉCESSAIRES 
# calcola_fst_pairwise_nobootstrap

################################ SLIDING WINDOW ################################

maf <- 0

#### Préparation des bornes à chaque itération : -------------------------------

slide <- 250000
window <- 1000000

vec_pos <- getPOS(VCF1)
n_pop <- length(list_pop)
low_bound <- seq(1, max(vec_pos) + slide, slide)
upper_bound <- seq(window, max(vec_pos) + slide + window, slide)


#### Pour chaque fenêtre : -----------------------------------------------------

fst_pairwise_sliding <- c()     

for (i in 1:length(low_bound)) {
  fenetre = which(vec_pos > low_bound[i] & vec_pos < upper_bound[i])
  
# ------- #### ... s'il y a des SNPs dans la fenêtre de positions ...  ----------
  if (length(fenetre) > 0) {
    data_temp <- VCF1[fenetre,]         # On sélectionne cette fenêtre de SNPs
    
          #### On calcule le pairwise FST de Reynolds ... :
    fst_temp <- calcola_fst_pairwise_nobootstrap(data_temp, list_pop, maf) 
    fst_pairwise_sliding <- rbind(fst_pairwise_sliding, fst_temp)
  } 
# ------- #### Sinon on met des NaN : ------------------------------------------
  else {
    fst_pairwise_sliding <- rbind(fst_pairwise_sliding, rep(NaN, n_pop*(n_pop-1)/2))
  }
}

# ------------------------------------------------------------------------------
############################### POUR LES FST ###################################
# ------------------------------------------------------------------------------

# Mise en forme du jeu de données : 
positions_plot <- (upper_bound+low_bound)/2
finale_fst <- na.omit(cbind(positions_plot, fst_pairwise_sliding))
finale_fst <- as_tibble(finale_fst)
write_csv(finale_fst, "/shared/projects/multiwhaling/multiwhaling/plot/SFS/fst_pairwise_sliding.csv")


list_sliding <- list()

# Graphique comprenant tous les FST des populations deux à deux :
pairs <- combn(seq_along(list_pop), 2) 
for (k in seq_len(ncol(pairs))) {
  i <- pairs[1, k]
  j <- pairs[2, k]
  
  y_values <- finale_fst[[k + 1]]   # <- on extrait la colonne AVANT ggplot
  colnames(finale_fst)[k + 1] <- paste(i, j, sep = "_")
  
  plot_sliding <- tibble(
    positions_plot = finale_fst$positions_plot,
    fst = y_values) |>
    ggplot(aes(x = positions_plot, y = fst)) +
    geom_point(shape = 1, alpha = 0.6) + 
    geom_hline(yintercept = quantile(y_values, probs = 0.95),
               color = "limegreen", linetype = 2) +
    geom_hline(yintercept = quantile(y_values, probs = 0.99),
               color = "red", linetype = 2) +
    scale_y_continuous(limits = c(0,0.25), breaks = seq(0, 0.25, 0.1)) +
    theme_minimal() +
    labs(title = paste(names(list_pop[i]), names(list_pop[j]), sep = "_"),
         x = "Position",
         y = "Pairwise FST") + 
    theme(axis.title = element_text(size = 30,face = "bold"), 
          title = element_text(size = 10, face = "bold"))
  
    list_sliding[[k]] <- plot_sliding 
}

write_csv(finale_fst, "/shared/projects/multiwhaling/multiwhaling/plot/SFS/fst_pairwise_sliding.csv")
saveRDS(list_sliding, "/shared/projects/multiwhaling/multiwhaling/plot/SFS/list_sliding_fst.RDS")



# pairwise_SNP_FST <- wrap_plots(list_sliding) + plot_layout(axis_titles = "collect", 
#                                                            axes = "keep",
#                                                            guides = "keep",
#                                                         design = "12345
#                                                                   #6789
#                                                                   ##ABC
#                                                                   ###DE
#                                                                   ####F") 
# pdf("~/Desktop/sliding_FST.pdf", width = 20, height = 12)
# pairwise_SNP_FST
# dev.off()


