################### 23/03/26 ######################
############## Sliding window FST #################
############## Achille RIVOLTELLA #################

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
#         - VCF avec SNPs                                                                             #
#         - objet (liste) contenant les noms des pops et les individus associés                       #
#         - n° chromosome                                                                             #
#                                                                                                     #
# SORTIES de la fonction :                                                                            #
#         - Matrice des FST par paires par fenêtre (positions)                                        #
#         - Objet RDS contenant tous les manhattan plots 2 à 2                                        #
#         - Graph ploté des manhattan plots                                                           # 
#######################################################################################################


#### ARGUMENTS :  
args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]
#chrom <- 21

window <- args[2]
slide <- args[3]
#window <- 500000
#slide <- 100000

VCF1 = read.vcfR(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/VCF_filtered_", 
                      chrom, ".vcf.gz", sep = ""))   
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/list_pop.RDS")


fst_sliding_window <- function(VCF, list_pop, window, slide){
  
  # Préparation des fenêtres : ---------------------------------------------------------
  vec_pos <- getPOS(VCF)
  low_bound <- seq(1, max(vec_pos) + slide, slide)
  upper_bound <- seq(window, max(vec_pos) + slide + window, slide)
  positions_plot <- (upper_bound + low_bound)/2
  n_pop <- length(list_pop)
  
  # ------------------------------------------------------------------------------------
  ################################## Calcul des FST ####################################
  # ------------------------------------------------------------------------------------
  fst_pairwise_sliding <- c()     
  
  pairs <- combn(seq_along(list_pop), 2)  # Création des paires de pop pour lesquelles on calcule le FST
  
  for (i in 1:length(low_bound)) {
    fenetre = which(vec_pos > low_bound[i] & vec_pos < upper_bound[i])
    
    # ------- #### ... s'il y a des SNPs dans la fenêtre de positions ...  ----------
    if (length(fenetre) > 0) {
      data_temp <- VCF1[fenetre,]         # On sélectionne cette fenêtre de SNPs
      
      #### On calcule le pairwise FST de Hudson ... :
      fst_temp <- data.frame(matrix(NA, nrow = 1, ncol = ncol(pairs)))
      
      for (k in seq_len(ncol(pairs))) {
        i <- pairs[1, k]
        j <- pairs[2, k]
        
        list_pop_hudson <- list(list_pop[[i]], list_pop[[j]])
        fst_pairwise <- fst.hudson(data_temp, list_pop_hudson)
        
        fst_temp[k] <- unlist(fst_pairwise[1])                      # Valeur moyenne de fst entre 2 pop
      }
      
      fst_pairwise_sliding <- rbind(fst_pairwise_sliding, fst_temp)
    } 
    # ------- #### Sinon on met des NaN : ------------------------------------------
    else {
      fst_pairwise_sliding <- rbind(fst_pairwise_sliding, rep(NaN, n_pop*(n_pop-1)/2))
    }
  }
  
  # Mise en forme du jeu de données : 
  finale_fst <- na.omit(cbind(positions_plot, fst_pairwise_sliding))
  finale_fst <- as_tibble(finale_fst)
  write_csv(finale_fst, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/FST/", chrom, "_fst_pairwise_sliding_hudson.csv", sep = ""))
  
  
  #### MISE EN FORME ET GRAPHIQUES : ---------------------------------------------------------------
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
  
  saveRDS(list_sliding, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/FST/", chrom, "_list_sliding_fst_hudson.RDS", sep = ""))
  
  
  pdf(paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/07_scan_genom/FST/", chrom, "_sliding_FST_hudson.pdf", sep = ""), 
      width = 20, height = 12)
  wrap_plots(list_sliding) + plot_layout(axis_titles = "collect", 
                                         axes = "keep",
                                         guides = "keep",
                                         design = "12345
                                                   #6789
                                                   ##ABC
                                                   ###DE
                                                   ####F")
  dev.off()
  
}
#########################################################################################################

fst_sliding_window(VCF1, list_pop, window, slide)






