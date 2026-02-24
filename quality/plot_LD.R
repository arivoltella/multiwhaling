#################### 24/02/26 #####################
##########  Plotting LD on a whole CHR ############
############# Achille RIVOLTELLA ##################

# Libraries : 
library(vcfR)
library(LDheatmap)

# Import data : 
VCF1 <- read.vcfR("/shared/projects/multiwhaling/Achille/VCF/sample_SNP.vcf.gz")
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/list_pop.RDS")

# Enlever les 2 individus apparentés : 
VCF1 <- VCF1[,c("FORMAT", unlist(list_pop))]
pos <- getPOS(VCF1)

# Et on ne garde qu'un SNP tous les 100 pour diminuer temps de calcul : 
pos_matrix <- seq(1,length(pos),100)

VCF_to_LDheatmap <- VCF1[pos_matrix,]


#### Visualiser le déséquilibre de liaison sur une heatmap : 
snpMatrix <- vcfR2SnpMatrix(VCF_to_LDheatmap)
pos_LDhm <- snpMatrix$genetic.distances

# On modifie le nom des colonnes pour que les SNP aient des noms : 
colnames(snpMatrix$data) <- paste0("SNP_", seq_len(ncol(snpMatrix$data)))

print("Type de l'objet que l'on donne à LD_heatmap")
str(snpMatrix$data)

# Sauvegarde du plot : 
png("LD_chr9.png")
LDheatmap(snpMatrix$data, 
          genetic.distances = pos_LDhm, 
          color = colorRampPalette(c("darkred","red", "white"))(20), 
          flip = T)
dev.off()

