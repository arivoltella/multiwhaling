#################### 19/03/26 #####################
############## Filtering of SNP ###################
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(ggplot2)
library(tidyverse)



#######################################################################################################
# ARGUMENTS à mettre en entrée :                                                                      #
#         - VCF avec SNPs ; VCF avec toutes les positions ;                                           #
#         - n° chromosome                                                                             #
#                                                                                                     #  
#                                                                                                     #
# SORTIES de la fonction :                                                                            #
#         - VCF filtré sur profondeur, Hz, NA, positions non polymorphiques                           #
#         - vecteur positions & VCF_all : mêmes filtres mais on garde les sites non polymorphiques    #
#         - names_ind et list_pop : association individus et pop (sous forme de liste ou tibble)      #
#         - VCF filtré sur la maf (0.05)                                                              #
#######################################################################################################

#### ARGUMENTs : 
args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]
#chrom <- 21

VCF <- read.vcfR(paste("/shared/projects/multiwhaling/VCF_AllPops/HumpbackAllIndiv_", chrom, "_GATK_TAG_Flowqual_Noindels_Norepeat_SNP.vcf.gz", sep = ""))
VCF_all <- read.vcfR(paste("/shared/projects/multiwhaling/VCF_AllPops/HumpbackAllIndiv_", chrom, "_GATK_TAG_Flowqual_Noindels_Norepeat.vcf.gz", sep = ""))


filtering_VCF <- function(VCF) {
  DP1 <- extract.gt(VCF, element='DP', as.numeric = TRUE) 
  name_object <- deparse(substitute(VCF))
  
  # -------------------------------------------------------------------------------------
  ################################ NOMS & POPULATIONS ###################################
  # -------------------------------------------------------------------------------------
  order <- c("BER", "KAR", "CHI", "PER", "WPA", "MAD", "SPM", "GOM")
  names_ind <- DP1 |> colnames() |> as_tibble()
  names_ind <- names_ind |>
    mutate(Population = str_extract(colnames(DP1), "^[A-Za-z]+|^[0-9]+")) |>
    mutate(Population = ifelse(Population == "A", "SPM",
                        ifelse(Population == "Ber", "BER", 
                        ifelse(Population == "Chi", "CHI",
                        ifelse(Population == "Kar", "KAR", 
                        ifelse(Population == "Pe", "PER",
                        ifelse(Population == "AN", "WPA", 
                        ifelse(Population == "MnD", "GOM",
                        ifelse(Population == "MnS", "GOM", "MAD"))))))))) |>
    filter(!(value %in% c("Chi2", "A20_13"))) |>
    rename(Individu = value) |>
    mutate(Population = factor(Population, levels = order)) |>
    arrange(Population)
  ############ Filtrer tous les individus détectés mauvais ##############
  list_pop <- split(names_ind$Individu, names_ind$Population)
  
  saveRDS(list_pop, "/shared/projects/multiwhaling/multiwhaling/whole_genome/data/list_pop.RDS")
  saveRDS(names_ind, "/shared/projects/multiwhaling/multiwhaling/whole_genome/data/names_ind.RDS")
  
  # Enlever les individus problématiques (relatedness) : 
  VCF <- VCF[,c("FORMAT", unlist(list_pop))] 
  
  # ----------------------------------------------------------------------------------
  ################################## FILTRES : #######################################
  # ----------------------------------------------------------------------------------
  ############ Choisir les filtres en fonction des plots de qualité obtenus ##########
  ################# Et des analyses que l'on veut faire ensuite ######################
  
  #### PROFONDEUR : -----------------------------------------------------------------------
  
  DP1 <- extract.gt(VCF, element='DP', as.numeric = TRUE) 
  depth = data.frame(depth_pos = apply(DP1, 1, mean,na.rm=T))
  
  VCF_DP <- subset(VCF, depth$depth_pos > 10 & depth$depth_pos < 60)        # À modifier en f° des besoins
  # Position avec trop ou pas assez de profondeur filtrées                   
  print("Depth of coverage filtered")
  
  #### HÉTÉROZZYGOTIE : -------------------------------------------------------------------
  
  geno1 <- as.data.frame(extract.gt(VCF_DP, element="GT", mask=F,as.numeric=F,
                                    return.alleles = F,convertNA = F,extract = T))
  # Génotype à chaque position pour chaque ind.
  n_ind <- dim(geno1)[2]
  het <- data.frame(het_pos = rowSums(geno1 == "0/1"))    # Nb d'Hz / position
  VCF_DP_hz <- subset(VCF_DP, het$het_pos < (n_ind*8)/10) # À modifier en f° des besoins
  print("Heterozigosity filtered")
  
  if (name_object == "VCF") {
    #### FOR VCF WITH ONLY SNPs
    #### ENLEVER POSITIONS MONOMORPHIQUES : --------------------------------------------------
    
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
    
    positions <- getPOS(VCF_DP_hz)
    snps_tot <- length(positions)
    genotypes <- extract.gt(VCF_DP_hz, element = "GT", mask = FALSE, as.numeric=F,
                            return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, 
                            convertNA = FALSE) 
    
    # Génotype pour chaque indiv à chaque pop : 
    count_allele <- apply(genotypes, 2, fun_geno_allele) |>
      as.numeric() |> matrix(ncol = n_ind, nrow = snps_tot, byrow=F)
    
    # À chaque position : 
    freq_ALT <- apply(count_allele, 1, sum, na.rm = T)/ (n_ind*2)   ## Fréquence de l'allèle alternatif 
    freq_REF <- 1 - freq_ALT                                        ## Fréquence de l'allèle de Réf
    
      ################################### POUR VCF NORMAL : ####################################
      # Enlève les sites non polymorphiques : 
      VCF_DP_hz_SNP <- subset(VCF_DP_hz, freq_REF > 0 | freq_REF < 1) 
    
      #### MISSING DATA : ----------------------------------------------------------------------
      genotypes <- extract.gt(VCF_DP_hz_SNP, element = "GT", mask = FALSE, 
                            as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, 
                            extract = TRUE, convertNA = FALSE)
      positions <- getPOS(VCF_DP_hz_SNP)
    
      NAs_pos <- rowSums(genotypes == "./.")
      VCF_DP_hz_SNP_NApos <- subset(VCF_DP_hz_SNP, NAs_pos < n_ind*0.2)
    
    
      #### SAUVEGARDE : ------------------------------------------------------------------------
      # Assigner chaque individu à une pop et les mettre dans le bon ordre 
      VCF_DP_hz_SNP_NA_ordered <- VCF_DP_hz_SNP_NApos[,c("FORMAT", unlist(list_pop))]       
    
      vcfR::write.vcf(VCF_DP_hz_SNP_NA_ordered, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/", 
                                                    name_object, "_filtered_", chrom, ".vcf.gz", sep = ""))
      print("Saved VCF_filtered")
      
      ############################## POUR VCF FILTRÉ SUR LA MAF : ###############################
      # On reprend le VCF déjà filtré sur DP et Hz : (on reprend ) 
      freq_both <- cbind(freq_ALT, freq_REF)
      maf <- apply(freq_both, 1, min)
    
      # On enlève aussi les sites dont la MAF < 0.05 : ----------------------------------------
      VCF_DP_hz_SNP_MAF <- subset(VCF_DP_hz_SNP, maf > 0.05) 
    
    
      #### MISSING DATA : ---------------------------------------------------------------------
      genotypes <- extract.gt(VCF_DP_hz_SNP_MAF, element = "GT", mask = FALSE, 
                            as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, 
                            extract = TRUE, convertNA = FALSE)
      positions <- getPOS(VCF_DP_hz_SNP_MAF)
    
      NAs_pos <- rowSums(genotypes == "./.")
      VCF_DP_hz_SNP_MAF_NApos <- subset(VCF_DP_hz_SNP_MAF, NAs_pos < n_ind*0.2)
    
    
      #### SAUVEGARDE : -----------------------------------------------------------------------
      # Assigner chaque individu à une pop et les mettre dans le bon ordre 
      VCF_DP_hz_SNP_MAF_NA_ordered <- VCF_DP_hz_SNP_MAF_NApos[,c("FORMAT", unlist(list_pop))]       
    
      vcfR::write.vcf(VCF_DP_hz_SNP_MAF_NA_ordered, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/", 
                                                    name_object, "_filtered_maf_", chrom, ".vcf.gz", sep = ""))
      print("Saved VCF_filtered_maf")
      
  } else
  { #### FOR VCF with ALL POSITIONS : 
    #### MISSING DATA : ---------------------------------------------------------------------
    genotypes <- extract.gt(VCF_DP_hz, element = "GT", mask = FALSE, 
                            as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, 
                            extract = TRUE, convertNA = FALSE)
    positions <- getPOS(VCF_DP_hz)
    
    NAs_pos <- rowSums(genotypes == "./.")
    VCF_DP_hz_SNP_NApos <- subset(VCF_DP_hz, NAs_pos < n_ind*0.2)
    
    #### SAUVEGARDE : -----------------------------------------------------------------------
    # Assigner chaque individu à une pop et les mettre dans le bon ordre 
    VCF_DP_hz_SNP_NA_ordered <- VCF_DP_hz_SNP_NApos[,c("FORMAT", unlist(list_pop))] 
    
    vec_pos <- getPOS(VCF_DP_hz_SNP_NA_ordered)
    saveRDS(vec_pos, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/vec_pos/vec_pos_all_", 
                           chrom, ".RDS", sep = ""))
    vcfR::write.vcf(VCF_DP_hz_SNP_NA_ordered, paste("/shared/projects/multiwhaling/multiwhaling/whole_genome/data/chrom/", 
                                                    name_object, "_filtered_", chrom, ".vcf.gz", sep = ""))
    print("Saved VCF_all_filtered")
  }
}
################################################################################################################

filtering_VCF(VCF)
filtering_VCF(VCF_all)


