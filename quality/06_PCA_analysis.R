################### 27/01/26 ######################
################# PCA analysis ####################
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)
library(pegas)
library(LEA)


#### Read the data #  
VCF1 = read.vcfR("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_maf_LD_Pruned.vcf.gz")
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/list_pop.RDS")

names_ind <- enframe(list_pop, name = "Population", value = "Individu") |>
  unnest(Individu)

#------------------------------------------------------------------------------------------------------
########## PCA et cluster analysis ########## ---------------------------------------------------------
#------------------------------------------------------------------------------------------------------

#### Transformer le fichier VCF en .GENO
fun_geno_allele<-function(data){
  for (i in 1:length(data)) {
    if (data[i]=="./.") {data[i]<-"9"}
    if (data[i]=="0/0") {data[i]<-0}
    if (data[i]=="0/1") {data[i]<-1}
    if (data[i]=="1/1") {data[i]<-2}
  }
  return(data)
}

genotype <- extract.gt(VCF1, element = "GT", mask = FALSE, as.numeric=F, 
                       return.alleles = FALSE, IDtoRowNames = TRUE, 
                       extract = TRUE, convertNA = FALSE) 

n_ind <- ncol(genotype)
snps_tot <- nrow(genotype)
qqq <- apply(genotype, 1, fun_geno_allele)
genotype_num <- matrix(as.numeric(qqq), ncol = n_ind, nrow = snps_tot, byrow = T)

write.table(genotype_num, file = "/shared/projects/multiwhaling/multiwhaling/data/VCF_maf_LD_Pruned.geno", quote = F, 
            sep = "", row.names = F, col.names=F)

#### Faire la PCA
whale_pca <- pca(input.file = "/shared/projects/multiwhaling/multiwhaling/data/VCF_maf_LD_Pruned.geno")
summary(whale_pca)

# Eigenvalues : 
pdf("/shared/projects/multiwhaling/multiwhaling/plot/PCA/eigenvalues_all_pop_maf_LD_pruned.pdf")
plot(whale_pca, lwd = 5, col = "blue", cex = .7, xlab = ("Factors"), ylab = "Eigenvalues")
dev.off()

# On récupère les données dans un nouveau dossier généré automatiquement par la fonction 'pca' utilisée
eigen_pca <- read.table("/shared/projects/multiwhaling/multiwhaling/data/VCF_maf_LD_Pruned.pca/VCF_maf_LD_Pruned.eigenvalues", header = F)
eigen_pca <- (eigen_pca / ((sum(eigen_pca)))) * 100
eigen_pca <- unlist(round(eigen_pca, digits = 2))

# Récupération des coordonnées de chaque individu dans les projections : 
proj <- read.table("/shared/projects/multiwhaling/multiwhaling/data/VCF_maf_LD_Pruned.pca/VCF_maf_LD_Pruned.projections", header = F)
proj <- as_tibble(proj) |>
  dplyr::select(1:4) |>
  mutate(pop = names_ind |> arrange(Population))
saveRDS(proj, "/shared/projects/multiwhaling/multiwhaling/plot/PCA/projections_maf_LD_Pruned.RDS")

pdf(paste("/shared/projects/multiwhaling/multiwhaling/plot/PCA/PCA_all_pop_maf_LD_Pruned.pdf"))
proj |>
  ggplot(aes(x = V1, y = V2, color = pop$Population)) + 
  geom_point(size = 2) + 
  theme_bw() + 
  labs(x = paste("PC1 ", as.character(eigen_pca[1]), "%"), 
       y = paste("PC2 ", as.character(eigen_pca[2]), "%"), 
       color = "Population")

proj |>
  ggplot(aes(x = V2, y = V3, color = pop$Population)) + 
  geom_point(size = 2) + 
  theme_bw() + 
  labs(x = paste("PC2 ", as.character(eigen_pca[2]), "%"), 
       y = paste("PC3 ", as.character(eigen_pca[3]), "%"), 
       color = "Population")

proj |>
  ggplot(aes(x = V3, y = V4, color = pop$Population)) + 
  geom_point(size = 2) + 
  theme_bw() + 
  labs(x = paste("PC3 ", as.character(eigen_pca[3]), "%"), 
       y = paste("PC4 ", as.character(eigen_pca[4]), "%"), 
       color = "Population")

dev.off()
