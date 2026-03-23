#################### 19/03/26 #####################
###### Estimation of kinship btw individuals ######
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)
library(pheatmap)

#######################################################################################################
# ARGUMENTS à mettre en entrée :                                                                      #
#         - .king et .king.id issus de PLINK2                                                         #
#                                                                                                     #
# SORTIES de la fonction :                                                                            #
#         - .csv des indices de relatedness entre indivs                                              #
#         - Heatmap de la relatedness entre indivs                                                    # 
#######################################################################################################

# À FAIRE SUR UN SEUL CHROMOSOME (PAS BESOIN DE PLUS)

# DANS SCRIPT BASH : 
# Calcul de la relatedness directement sur le cluster avec PLINK2 : 
# plink2 --vcf VCF_de_base.vcf.gz --make-king --out data/relatedness/matrix_all

# -------------------------------------------------------------------------------------
################################# Relatedness : #######################################
# -------------------------------------------------------------------------------------

# Importer l'output de la kinship : 
king_id <- read.table("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/02_relatedness/matrix_all.king.id")
king <- readLines("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/02_relatedness/matrix_all.king")


relatedness <- function(king_id, king) {
  
  relatedness_king <- c()
  
  for (i in 1:length(king)) {
    
    # Récupération des valeurs individuelles et ajout de NA sur le reste de la ligne : 
    ind <- as.numeric(unlist(strsplit(king[i], split = "\t")))
    relatedness <- c(ind, rep(NA, length(king)- i))
    
    # Final matrix 
    relatedness_king <- rbind(relatedness_king, relatedness)
  }
  relatedness_king <- rbind(rep(NA, length(king)),relatedness_king)    # Ajouter une ligne de NA
  relatedness_king <- cbind(relatedness_king,rep(NA, length(king_id))) # Ajouter une colonne de NA
  
  colnames(relatedness_king) <- king_id$V1 
  rownames(relatedness_king) <- king_id$V1        # Ajout des noms d'individus 
  relatedness_king <- as.data.frame(relatedness_king)
  saveRDS(relatedness_king, "/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/02_relatedness/kinship.RDS")
  
  png("/shared/projects/multiwhaling/multiwhaling/whole_genome/plot/02_relatedness/heatmap_kinship.png", 
      width = 850, height = 750)
  pheatmap(relatedness_king, cluster_rows=F, cluster_cols=F, na_col="white",main = "Pairwise relatedness", 
                              color = colorRampPalette(c("seashell1", "yellow", "firebrick3"))(60), 
                              breaks = seq(from = -0.1, to = 0.5, by = 0.01))
  dev.off()
}
##########################################################################################################


relatedness(king_id, king)