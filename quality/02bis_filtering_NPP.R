#################### 27/01/26 #####################
############## Filtering of NPP ###################
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)

#### Read the VCF #  
VCF_all = read.vcfR("/shared/projects/multiwhaling/Achille/VCF/HumpbackTot_9_GATK_TAG_Flowqual_Noindels_Norepeat.vcf.gz")    # VCF avec toutes les positions
DP1 <- extract.gt(VCF_all, element='DP', as.numeric = TRUE) 

list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/list_pop.RDS")


# ---------------------------------------------------------------------------------------
# VCF toutes positions : SNP & NPP ---------------------------------------------------
# ---------------------------------------------------------------------------------------


#### PROFONDEUR : ----
depth = data.frame(depth_pos = apply(DP1, 1, mean, na.rm=T))

VCF_DP <- subset(VCF_all, depth$depth_pos > 10 & depth$depth_pos < 60)        # À modifier en f° des besoins
# Positions avec trop ou pas assez de profondeur filtrées
rm(VCF_all, depth)

#### HÉTÉROZZYGOTIE : ---- 
geno1 <- as.data.frame(extract.gt(VCF_DP, element="GT", mask=F,as.numeric=F,return.alleles = F,
                                  convertNA = F,extract = T))
# Génotype à chaque position pour chaque ind.
n_ind <- dim(geno1)[2]

het <- data.frame(het_pos = rowSums(geno1 == "0/1"))      # Nb d'Hz / position

VCF_DP_hz <- subset(VCF_DP, het$het_pos < (n_ind*8)/10)      #### À modifier en f° des besoins
rm(VCF_DP, geno1, het)          ## On supprime les objets intermédiaires pour éviter de prendre trop de mémoire pdt le calcul

# #### MAF : on fait pas car on veut garder TOUS les sites 


#### NA : ----
genotypes <- extract.gt(VCF_DP_hz, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
positions <- getPOS(VCF_DP_hz)

NAs_pos <- rowSums(genotypes == "./.")
VCF_DP_hz_NA <- subset(VCF_DP_hz, NAs_pos < n_ind*0.2)
rm(VCF_DP_hz, genotypes, positions, NAs_pos)

VCF_DP_hz_SNP_NA_ordered <- VCF_DP_hz_NA[,c("FORMAT", unlist(list_pop))]       
# Assigner chaque individu à une pop et les mettre dans le bon ordre 
rm(VCF_DP_hz_NA)

saveRDS(VCF_DP_hz_SNP_NA_ordered, "/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered_all.RDS")

