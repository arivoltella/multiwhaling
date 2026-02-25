#################### 24/02/26 #####################
##########  Plotting LD on a whole CHR ############
############# Achille RIVOLTELLA ##################

# Libraries : 
library(vcfR)
library(LDheatmap)

# Import data : 
VCF1 <- read.vcfR("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_chr9.vcf.gz")


# Et on ne garde qu'un SNP tous les 100 pour diminuer temps de calcul :
pos <- getPOS(VCF1)
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
png("plot/LD_heatmap/LD_chr9.png")
LDheatmap(snpMatrix$data, 
          genetic.distances = pos_LDhm, 
          color = colorRampPalette(c("darkred","red", "white"))(20), 
          flip = T)
dev.off()

