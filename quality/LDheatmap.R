#################### 24/02/26 #####################
##########  Plotting LD on a whole CHR ############
############# Achille RIVOLTELLA ##################

# Version Hugo et Julien 

# Import libraries : 
library(vcfR)
library(LDheatmap)
library(RColorBrewer)
library(snpStats)

#### Functions  ####
#Bin function = extrait un SNP tout les certains nombres définis et fait un vcf avec
bin_snps<-function(vcf,bin_size){
  res <- c()
  pos <- getPOS(vcf)
  seq_I <- seq(1,length(pos),bin_size)
  for (i in seq_I){
    binf <- i
    bsup <- i+bin_size-1
    pos_bin <- sample(binf:bsup,1)
    
    if(pos_bin > length(pos)) { break }
    
    res<-c(res,pos_bin)
  }
  print(res)
  vcf_clean <- vcf[res,]
  return(vcf_clean)
}

#fonction extract matrix coordinates
get_matrix_coordinates <- function(index, nrow) {
  row <- (index - 1) %% nrow + 1
  col <- (index -1) %/% nrow + 1
  return(c(row, col))
}

#extract LD function = extrait les cordonnées des SNP au dessus d'un certain seuil de corrélation (R2)
extract_LD <- function(LD_objet, seuil){
  res = c()
  s = which(LD_objet$LDmatrix >= seuil)
  for (i in 1:length(s)){
    tt = get_matrix_coordinates(index = s[i], nrow = nrow(LD_objet$LDmatrix))
    pos1 = LD_objet$genetic.distances[tt[1]]
    pos2 = LD_objet$genetic.distances[tt[2]]
    res = rbind(res, c(pos1,pos2))
    colnames(res) = c("snp1", "snp2")
  }
  return(res)
}




################## ANALYSIS #####################
vcf = "/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_chr9.vcf.gz"
vcf1<-read.vcfR(vcf, verbose = TRUE)
vcfname <- gsub(".vcf.gz", "", basename(vcf))

poubelle_taille = 100

##### Boucle pour faire les vcf par chromosome : #### --------------------------
chromosomes <- unique(vcf1@fix[, "CHROM"])

for (chr in chromosomes) {
  
  vcf_sub <- vcf1[vcf1@fix[, "CHROM"] == chr, ]
  
  nom_fichier <- paste0("chr_",chr, ".vcf")

  assign(nom_fichier, vcf_sub)
}

rm(vcf1,vcf_sub)

##### Réalisation des plots de LD : ##### --------------------------------------
for (chr in chromosomes) {
# 1/ bin data
vcf_file <- paste0("chr_",chr, ".vcf")  #ajouter l'extension
data_bin <- bin_snps(vcf = get(vcf_file), bin_size = poubelle_taille)
print(data_bin) #regarder s'il y a bien 100 fois moins de SNP

# 2/ run LDheatmap
snp_matrix <- vcfR2SnpMatrix(data_bin, phased = NULL, subjects = NULL)

# Nom du fichier png et paramètre du png
nompng = paste0("LDheatmap_",vcfname,"_",chr,".png")
png(nompng, width = 11.69*1000, height = 8.27*1000, res = 1000)

info_titre=paste0("fichier = ",vcfname," / bin size = ",poubelle_taille)
LD <- LDheatmap(snp_matrix$data, snp_matrix$genetic.distances,LDmeasure = "r", flip = T, 
                color = colorRampPalette(c("darkred", "lightcoral", "white"))(50),title = info_titre)

dev.off()

##### 3/ extract coordinates of LD snps ##### ----------------------------------
coord = extract_LD(LD_objet = LD, seuil = 0.9)
write.table(coord, paste0("plot/LD_heatmap/Coord_",vcfname,"_",chr,".txt"), sep = " ", row.names = TRUE, col.names = TRUE)
}

