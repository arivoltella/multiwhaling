#################### 24/02/26 #####################
##########  Plotting LD on a whole CHR ############
############# Achille RIVOLTELLA ##################

# Libraries : 
library(vcfR)
library(LDheatmap)

# Import data : 
VCF1 <- read.vcfR("/shared/projects/multiwhaling/Achille/VCF/HumpbackTot_9_GATK_TAG_Flowqual_Noindels_Norepeat_SNP.vcf.gz")
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/list_pop.RDS")

# Enlever les 2 individus apparentés : 
VCF1 <- VCF1[,c("FORMAT", unlist(list_pop))]   

# Visualiser le déséquilibre de liaison sur une heatmap : 
snpMatrix <- vcfR2SnpMatrix(VCF1)
pos <- snpMatrix$genetic.distances

png("/shared/projects/multiwhaling/multiwhaling/plot/LD_heatmap/LD_chr9.png")
LDheatmap(snpMatrix$data, 
          genetic.distances = pos, 
          color = colorRampPalette(c("red","white"))(20), 
          flip = T)
dev.off()

