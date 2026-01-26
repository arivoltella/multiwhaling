#################### 23/01/26 #####################
##### Population structure of Humpback whales #####
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)
library(pegas)
library(coala)
library(MASS)
library(LEA)

#### Read the VCF #### 
VCF1 = read.vcfR("../Achille/data/sample_SNP.vcf.gz")
DP1 <- extract.gt(VCF1, element='DP', as.numeric = TRUE) 


################# FILTRES : ################ ----------------------------------------------

#### PROFONDEUR : 
depth = data.frame(depth_pos = apply(DP1, 1, mean,na.rm=T))

VCF_DP <- subset(VCF1, depth$depth_pos > 10 & depth$depth_pos < 50)
# Position avec trop ou pas assez de profondeur filtrées

#### HÉTÉROZZYGOTIE : 
geno1 <- as.data.frame(extract.gt(VCF_DP, element="GT", mask=F,as.numeric=F,return.alleles = F,
                                  convertNA = F,extract = T))
# Génotype à chaque position pour chaque ind.
n_ind <- dim(geno1)[2]

het <- data.frame(het_pos = rowSums(geno1 == "0/1"))# Nb d'Hz / position

VCF_DP_hz <- subset(VCF_DP, het$het_pos < (n_ind*8)/10)

#### MAF :
# Fonction qui assigne les génotypes : 
fun_geno_allele<-function(data){
  for (i in 1:length(data)) {
    if (data[i]=="./.") {data[i]<-"NA"}
    if (data[i]=="0/0") {data[i]<-0}
    if (data[i]=="0/1") {data[i]<-1}
    if (data[i]=="1/1") {data[i]<-2}
  }
  return(data)
}
# assigne hétéro/homozygotie : 
fun_geno_mod<-function(data){
  for (i in 1:length(data)) {
    if (data[i]=="./.") {data[i]<-"NA"}
    if (data[i]=="0/0") {data[i]<-0}
    if (data[i]=="0/1") {data[i]<-1}
    if (data[i]=="1/1") {data[i]<-0}
  }
  return(data)
}

positions <- getPOS(VCF_DP_hz)
snps_tot <- length(positions)

genotypes <- extract.gt(VCF_DP_hz, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0

# Génotype pour chaque indiv à chaque pop : 
count_allele <- apply(genotypes, 2, fun_geno_allele) |>
  as.numeric() |> matrix(ncol = n_ind, nrow = snps_tot, byrow=F)

# À chaque position : 
freq_ALT <- apply(count_allele, 1, sum)/ (n_ind*2)      ## Fréquence de l'allèle alternatif 
freq_REF <- 1 - freq_ALT                                  ## Fréquence de l'allèle de Réf
freq_both <- cbind(freq_ALT, freq_REF)

VCF_DP_hz_SNP <- subset(VCF_DP_hz, freq_REF > 0 | freq_REF < 1)
# Enlève les sites non polymorphiques 
saveRDS(VCF_DP_hz_SNP, "data/VCF_filtered.RDS")


#### NOMS & POPULATIONS #### ---------------------------------------------------
names <- colnames(VCF1@gt)[-1]
namest <- as_tibble(names)
names_ind <- namest |>
  mutate(pop = str_extract(names, "^[A-Za-z]+|^[0-9]+")) |>
  mutate(pop = ifelse(pop == "A", "N_ATL",
                      ifelse(pop == "Ber", "BERING", 
                             ifelse(pop == "Chi", "CHILI", 
                                    ifelse(pop == "Kar", "KARAGINSKY", 
                                           ifelse(pop == "Pe", "PEROU", "MADAGASCAR"))))))
list_pop <- split(names_ind$value, names_ind$pop)
saveRDS(list_pop, "data/list_pop.RDS")


##### FIS ------------------------------------------------------------------

genotypes <- extract.gt(VCF1, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
positions <- getPOS(VCF1)
snps_tot <- length(positions)

count_het <- apply(genotypes, 2, fun_geno_mod) |>
  as.numeric() |> matrix(ncol=n_ind, nrow=snps_tot, byrow=F) 
### matrice numerique 0/1

count_allele <- apply(genotypes, 2, fun_geno_allele) |>
  as.numeric() |> matrix(ncol = n_ind, nrow = snps_tot, byrow=F)
### matrice numerique 0/1/2

# Boucle permettant de calculer Hz attendue et observée : 
a=0
vec_het_obs <-c()
vec_het_exp <-c()

for (i in seq_along(list_pop)) {
  
  ncols <- length(list_pop[[i]])
  col_range <- (a+1):(a+ncols)
  
  hetero_obs <- apply(count_het[,col_range], 1, sum)/ncols
  hetero_exp_temp <- apply(count_allele[,col_range], 1, sum)/(ncols*2)
  hetero_exp <- 2 * hetero_exp_temp * (1 - hetero_exp_temp)
  
  a <- a + ncols
  
  vec_het_obs <- cbind(vec_het_obs,hetero_obs)
  vec_het_exp <- cbind(vec_het_exp,hetero_exp)
}
vec_het_obs <- cbind(positions,vec_het_obs) # vector observed HET by pop
colnames(vec_het_obs) <- c("Position", names(list_pop))
vec_het_exp <- cbind(positions,vec_het_exp) # Vector expected HET by pop
colnames(vec_het_exp) <- c("Position", names(list_pop))

all_FIS<- c()       
for (i in seq_along(list_pop)) {
  all_FIS <- cbind(all_FIS,(vec_het_exp[,i+1] - vec_het_obs[,i+1])/(vec_het_exp[,i+1]))
}

colnames(all_FIS) <- c(names(list_pop))
all_FIS <- as_tibble(all_FIS) |>
 mutate(positions = positions) |>
  pivot_longer(cols = 1:6, 
               names_to = "pop", 
               values_to = "FIS")

all_FIS |>
  ggplot(aes(x = positions, y = FIS, color = pop))+
  geom_point()

all_FIS |>
  filter(!is.na(FIS)) |>
  reframe(mean_FIS = mean(FIS),
          .by = pop)



################ DIVERSITÉ GÉNÉTIQUE et SFS ---------------------------------------------------
# Div génétique et SFS pour tout le jeu de données 

### Fonctions sources 
source("./../Achille/functions_for_td.r")

data_bin_all <- vcfR2DNAbin(VCF_DP_hz_SNP, extract.indels = TRUE, consensus = FALSE,
  unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL, verbose = TRUE) 

sfs_tot_all <- site.spectrum(data_bin_all, folded = TRUE) #### sfs folded con funzione PEGAS
sfs_tot_all <- as.vector(sfs_tot_all)
barplot(sfs_tot_all) 

# Pour regarder si le SFS dévie du modèle neutre : 
norm_sfs_tot_all <- calcola_normalized_foldedSFS(sfs_tot_all) 
plot(norm_sfs_tot_all, type = "l", ylim = c(0, 0.5))

# ADD transformed spectrum of a Kingman pop #
# This is done to check the variation in the SFS class compared to Kingman's expectations
kingman_all <- coal_model(2*length(sfs_tot_all), loci_number = 30000, loci_length = 100, ploidy = 1) + 
  feat_mutation(1) + sumstat_sfs() 
sumstats <- simulate(kingman_all, seed = 23432)
sim_folded_sfs <- fold_dafSFS_coala(sumstats$sfs)
norm_simulated_sfs <- calcola_normalized_foldedSFS(sim_folded_sfs)
lines(norm_simulated_sfs, col = 4)

# Métriques de diversité génétique : 
div_all <- unlist(calcola_TD_folded(sfs_tot_all)) 
names(div_all) <- c("TD", "Pi", "W", "S")
div_all




#------------------------------------------------------------------------------------------------------
########## SFS et DIVERSITÉ GÉNÉTIQUE pop par pop ########## ------------------------------------------
#------------------------------------------------------------------------------------------------------

#### Data : 
list_pop <- readRDS("data/list_pop.RDS")
VCF_filtered <- readRDS("data/VCF_filtered.RDS")

#### Création des listes par population : 
for (i in seq_along(list_pop)) {
  
  # Tri des individus pour les assigner la pop à laquelle on s'intéresse
  pop <- list_pop[i]
  data <- VCF_filtered[,c("FORMAT", unlist(pop))] 
  
  # Calcul du SFS : 
  dati_bin <- vcfR2DNAbin(data, extract.indels = TRUE, consensus = FALSE, 
                          unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL,
                          verbose = TRUE)
  sfs_tot <- as.vector(site.spectrum(dati_bin, folded = TRUE)) 
  
  ######## Plot SFS  : 
  pdf(paste("plot/SFS/", names(list_pop[i]), "_sfs.pdf", sep = ""))
  barplot(sfs_tot) 
  dev.off()

  norm_sfs_tot <- calcola_normalized_foldedSFS(sfs_tot) 
  ## transformed SFS to inspect deviation from the neutral model
  ## ADD transformed spectrum of a Kingman pop ########### This is done to check the variation in the SFS class compared to Kingman's expectations
  kingman <- coal_model(2*length(sfs_tot), loci_number = 30000, loci_length = 100, ploidy = 1) + 
    feat_mutation(1) + sumstat_sfs() 
  sumstats <- simulate(kingman, seed = 23432)
  sim_folded_sfs <- fold_dafSFS_coala(sumstats$sfs)
  norm_simulated_sfs <- calcola_normalized_foldedSFS(sim_folded_sfs)
  
  ######### Plot normalized SFS + modèle neutre : 
  pdf(paste("plot/SFS/", names(list_pop[i]), "_norm_sfs.pdf", sep = ""))
  plot(norm_sfs_tot, type = "l", ylim = c(0, 0.5))
  lines(norm_simulated_sfs, col=4)
  dev.off()
  
  ######### Diversité génétique : Tajima's D, thetaPi, thetaS, S
  div_gen <- unlist(calcola_TD_folded(sfs_tot)) 
  names(div_gen) <- c("TD", "Pi", "W", "S")
  write.csv(div_gen, paste("plot/SFS/", names(list_pop[i]), "_div_gen.csv", sep = ""))
}


#------------------------------------------------------------------------------------------------------
########## FST des populations deux à deux ########## ------------------------------------------
#------------------------------------------------------------------------------------------------------

# Ré-ordonner les individus par pop : 
data <- VCF_filtered[,c("FORMAT", unlist(list_pop))]


fst_moy <- data.frame(matrix(NA, nrow = length(list_pop), ncol = length(list_pop), 
                             dimnames = list(names(list_pop), names(list_pop))))

pairs <- combn(seq_along(list_pop), 2)  # Création des paires de pop pour lesquelles on calcule le FST

for (k in seq_len(ncol(pairs))) {
  i <- pairs[1, k]
  j <- pairs[2, k]
  
  list_pop_hudson <- list(list_pop[[i]], list_pop[[j]])
  fst_pairwise <- fst.hudson(data, list_pop_hudson)
  
  fst_moy[i, j] <- fst_pairwise[1]                      # Valeur moyenne de fst entre 2 pop
  fst_moy[j, i] <- fst_pairwise[1]                      # pour matrice symétrique
    
  fst_per_SNP <- unlist(fst_pairwise[2])                # Fst par position SNP
    
  # Manhattan plot 
  manh <- cbind(rep(1, length(fst_per_SNP)), positions, fst_per_SNP)
  manh <- data.frame(names(fst_per_SNP), manh)
  colnames(manh) <- c("SNP","CHR","BP","P")
    
  pdf(paste("plot/FST/", names(list_pop[i]), names(list_pop[j]), "FST.pdf", sep = "_"))
  manhattan(manh, ylim = c(0,1), logp = FALSE, ylab = "FST", suggestiveline = F, genomewideline = F)
  dev.off()
  }



#------------------------------------------------------------------------------------------------------
########## PCA et cluster analysis ########## ---------------------------------------------------------
#------------------------------------------------------------------------------------------------------

# Data : 
data <- readRDS("data/VCF_filtered.RDS")

#### Transformer le fichier VCF en .GENO
fun_geno_allele<-function(data){
  for (i in 1:length(data)) {
    if (data[i]=="./.") {data[i]<-"NA"}
    if (data[i]=="0/0") {data[i]<-0}
    if (data[i]=="0/1") {data[i]<-1}
    if (data[i]=="1/1") {data[i]<-2}
  }
  return(data)
}

genotype <- extract.gt(data, element = "GT", mask = FALSE, as.numeric=F, 
                       return.alleles = FALSE, IDtoRowNames = TRUE, 
                       extract = TRUE, convertNA = FALSE) 

n_ind <- ncol(genotype)
snps_tot <- nrow(genotype)
qqq <- apply(genotype, 1, fun_geno_allele)
genotype_num <- matrix(as.numeric(qqq), ncol = n_ind, nrow = snps_tot, byrow = T)

write.table(genotype_num, file="data/VCF_filtered.geno", quote = F, 
            sep = "", row.names = F, col.names=F)

#### Faire la PCA
whale_pca <- pca(input.file = "data/VCF_filtered.geno")
summary(whale_pca)

# Eigenvalues : 
pdf(paste("plot/PCA/eigenvalues_all_pop.pdf"))
plot(whale_pca, lwd = 5, col = "blue", cex = .7, xlab = ("Factors"), ylab = "Eigenvalues")
dev.off()

# On récupère les données dans un nouveau dossier généré automatiquement par la fonction 'pca' utilisée
eigen_pca <- read.table("data/VCF_filtered.pca/VCF_filtered.eigenvalues", header = F)
eigen_pca <- (eigen_pca / ((sum(eigen_pca)))) * 100
eigen_pca <- unlist(round(eigen_pca, digits = 2))

# Récupération des coordonnées de chaque individu dans les projections : 
proj <- read.table("data/VCF_filtered.pca/VCF_filtered.projections", header = F)
proj <- as_tibble(proj) |>
  select(1:2) |>
  mutate(pop = names_ind |> arrange(pop))

pdf(paste("plot/PCA/PCA_all_pop.pdf"))
proj |>
  ggplot(aes(x = V1, y = V2, color = pop$pop)) + 
  geom_point(size = 3) + 
  theme_bw() + 
  labs(x = paste("PC1 ", as.character(eigen_pca[1]), "%"), 
       y = paste("PC2 ", as.character(eigen_pca[2]), "%"), 
                 color = "Population")
dev.off()


# ------------------------------------------------------------------
#### sNMF clustering  ----------------------------------------------
# ------------------------------------------------------------------

whale_snmf <- snmf(input.file= "data/VCF_filtered.geno", K= 1:6, repetitions = 3, project="new", entropy=T)
show(whale_snmf)
summary(whale_snmf)
plot(whale_snmf, col = "blue", pch = 19, cex = 1.2) ### Identifier le K avec l'entropie la plus faible
best <- which.min(cross.entropy(whale_snmf, K = 5))
my.colors <- c("royalblue", "seagreen4", "peru", "khaki", "indianred2")

adm_coeff <- Q(whale_snmf, K = 5, run = best) 
barplot(t(adm_coeff), names = names_ind$value, las = 2, cex.names = 0.6, border = NA, space = 0, 
        col = my.colors, 
        xlab = "Individuals", 
        ylab = "Ancestry proportions", 
        main = "Ancestry matrix")




