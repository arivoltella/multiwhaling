#################### 27/01/26 #####################
############## Filtering of SNP ###################
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)

#### Read the VCF #  
VCF1 <- read.vcfR("/shared/projects/multiwhaling/Achille/VCF/HumpbackTot_9_GATK_TAG_Flowqual_Noindels_Norepeat_SNP.vcf.gz")
DP1 <- extract.gt(VCF1, element='DP', as.numeric = TRUE) 
# /shared/projects/multiwhaling/Achille/VCF

# -------------------------------------------------------------------------------------
############## NOMS & POPULATIONS ################# -----------------------------------
# -------------------------------------------------------------------------------------

names <- colnames(VCF1@gt)[-1]
namest <- as_tibble(names)
names_ind <- namest |>
  mutate(pop = str_extract(names, "^[A-Za-z]+|^[0-9]+")) |>
  mutate(pop = ifelse(pop == "A", "N_ATL",
                      ifelse(pop == "Ber", "BERING", 
                             ifelse(pop == "Chi", "CHILI", 
                                    ifelse(pop == "Kar", "KARAGINSKY", 
                                           ifelse(pop == "Pe", "PEROU", "MADAGASCAR"))))))|>
  filter(!(value %in% c("Chi2", "A20_13")))
list_pop <- split(names_ind$value, names_ind$pop)
saveRDS(list_pop, "/shared/projects/multiwhaling/multiwhaling/data/list_pop.RDS")



# ----------------------------------------------------------------------------------
################# FILTRES : ################ ---------------------------------------
# ----------------------------------------------------------------------------------

############ Choisir les filtres en fonction des plots de qualité obtenus #########
################## Et des analyses que l'on veut faire ensuite ####################


#### PROFONDEUR : ----
depth = data.frame(depth_pos = apply(DP1, 1, mean,na.rm=T))
                                                                           
VCF_DP <- subset(VCF1, depth$depth_pos > 10 & depth$depth_pos < 60)        # À modifier en f° des besoins
# Position avec trop ou pas assez de profondeur filtrées                   


#### HÉTÉROZZYGOTIE : ---- 
geno1 <- as.data.frame(extract.gt(VCF_DP, element="GT", mask=F,as.numeric=F,return.alleles = F,
                                  convertNA = F,extract = T))
# Génotype à chaque position pour chaque ind.
n_ind <- dim(geno1)[2]

het <- data.frame(het_pos = rowSums(geno1 == "0/1"))      # Nb d'Hz / position

VCF_DP_hz <- subset(VCF_DP, het$het_pos < (n_ind*8)/10)               #### À modifier en f° des besoins


#### MAF : ----
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
freq_ALT <- apply(count_allele, 1, sum, na.rm = T)/ (n_ind*2)        ## Fréquence de l'allèle alternatif 
freq_REF <- 1 - freq_ALT                                  ## Fréquence de l'allèle de Réf
freq_both <- cbind(freq_ALT, freq_REF)
maf <- apply(freq_both, 1, min)

# Enlève les sites non polymorphiques : 
VCF_DP_hz_SNP <- subset(VCF_DP_hz, freq_REF > 0 | freq_REF < 1) 
VCF_DP_hz_SNP_MAF <- subset(VCF_DP_hz_SNP, maf > 0.05) 


#### NA : ----
genotypes <- extract.gt(VCF_DP_hz_SNP_MAF, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
positions <- getPOS(VCF_DP_hz_SNP_MAF)

NAs_pos <- rowSums(genotypes == "./.")
VCF_DP_hz_SNP_NApos <- subset(VCF_DP_hz_SNP_MAF, NAs_pos < n_ind*0.2)


VCF_DP_hz_SNP_NA_ordered <- VCF_DP_hz_SNP_NApos[,c("FORMAT", unlist(list_pop))]       
# Assigner chaque individu à une pop et les mettre dans le bon ordre 

saveRDS(VCF_DP_hz_SNP_NA_ordered, "/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_maf.RDS")
# Sauvegarder en .vcf aussi 
write.vcf(VCF_DP_hz_SNP_NA_ordered, "/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_maf.vcf.gz")
