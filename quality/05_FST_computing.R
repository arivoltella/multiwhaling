################### 27/01/26 ######################
############### Compute the FST ###################
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)
library(pegas)
library(pheatmap)



#### Read the data #  
VCF1 = read.vcfR("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_maf_LD_Pruned.vcf.gz")
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/list_pop.RDS")

source("/shared/projects/multiwhaling/multiwhaling/quality/functions_for_td.r")

#------------------------------------------------------------------------------------------------------
########## FST des populations deux à deux ########## -------------------------------------------------
#------------------------------------------------------------------------------------------------------

# Data nécessaires : 
positions <- getPOS(VCF1)

fst_moy <- data.frame(matrix(NA, nrow = length(list_pop), ncol = length(list_pop), 
                             dimnames = list(names(list_pop), names(list_pop))))

pairs <- combn(seq_along(list_pop), 2)  # Création des paires de pop pour lesquelles on calcule le FST

for (k in seq_len(ncol(pairs))) {
  i <- pairs[1, k]
  j <- pairs[2, k]
  
  list_pop_hudson <- list(list_pop[[i]], list_pop[[j]])
  fst_pairwise <- fst.hudson(VCF1, list_pop_hudson)
  
  fst_moy[i, j] <- unlist(fst_pairwise[1])                      # Valeur moyenne de fst entre 2 pop
  
  fst_per_SNP <- unlist(fst_pairwise[2])                # Fst par position SNP
  
  # Manhattan plot 
  manh <- cbind(rep(1, length(fst_per_SNP)), positions, fst_per_SNP)
  manh <- data.frame(names(fst_per_SNP), manh)
  colnames(manh) <- c("SNP","CHR","BP","P")

  pdf(paste("/shared/projects/multiwhaling/multiwhaling/plot/FST/", names(list_pop[i]), names(list_pop[j]), "FST_pruned.pdf", sep = "_"))
  manhattan(manh, ylim = c(0,1), logp = FALSE, ylab = "FST", suggestiveline = F, genomewideline = F)
  dev.off()
}

write.csv(fst_moy, file = "/shared/projects/multiwhaling/multiwhaling/plot/FST/mean_fst_pop_pruned.csv")

png("/shared/projects/multiwhaling/multiwhaling/plot/FST/fst_pairwise_heatmap_pruned.png", height = 750, width = 750)
pheatmap(fst_moy, cluster_rows=F, cluster_cols=F, na_col="white",main = "Pairwise Fst", 
         color = colorRampPalette(c("seashell1", "yellow", "firebrick3"))(50))
dev.off()