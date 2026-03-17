#################### 11/03/26 #####################
############ Sliding window SFS : #################
############# Achille RIVOLTELLA ##################

# Libraries : 
library(vcfR)
library(tidyverse)
library(rlist)
library(pegas)
library(ggplot2)
library(adegenet)
library(patchwork)

# Data : 
VCF1 <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered.RDS")
VectPosFiltered <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/vect_position/pos_all.RDS")         # Toutes positions du chromosomes 
list_pop <- readRDS("data/list_pop.RDS")

source("quality/functions_for_td.r")

# ATTENTION : IMPORTER LES FONCTIONS NÉCESSAIRES 
sliding_spectre_onepop<-function(dati, window, slide){
  dim_locus<-window
  interval<-slide
  vettore_snps_pos<-getPOS(dati)
  vettore_nomi<-c(names(dati@gt[1,])) #### creo vettore nomi dal vcf
  vettore_nomi<-vettore_nomi[-1]      #### creo vettore nomi dal vcf
  n_ind<-length(vettore_nomi)
  low_bound<-seq(1, max(vettore_snps_pos)+interval, interval)
  upper_bound<-seq(dim_locus, max(vettore_snps_pos)+interval+dim_locus, interval)
  per_plot<-(upper_bound+low_bound)/2
  
  ris_spectre_sliding<-c()
  for (i in 1:length(low_bound)) {
    sss=which(vettore_snps_pos > low_bound[i] & vettore_snps_pos < upper_bound[i] )
    if (length(sss)>0) {
      dati_temp<-dati[sss,]
      temp_spectre<-c()
      dati_bin<-vcfR2DNAbin(dati_temp, extract.indels = TRUE, consensus = FALSE,  unphased_as_NA = FALSE, ref.seq = NULL, start.pos = NULL, verbose = TRUE) ### trasformo in pratica il vcf in sequenze, e posso quindi calcolare sia il SFS che la nucleotide diversity. Ma non l'heterozygosity
      sfs_tot<-site.spectrum(dati_bin, folded = TRUE) #### sfs folded con funzione PEGAS
      temp_spectre<-c(temp_spectre,sfs_tot)
      ris_spectre_sliding<-rbind(ris_spectre_sliding,temp_spectre)
    } else {
      ris_spectre_sliding<-rbind(ris_spectre_sliding, rep(NaN,n_ind))
    }
  }
  
  finale_spectre<-cbind(per_plot,ris_spectre_sliding)
  finale_spectre<-na.omit(finale_spectre)
  
  return(finale_spectre)
}

#### FILTRER LES NA POUR ANALYSE SFS ####

################################ SLIDING WINDOW ################################
z = 1 # Bering

# Filtre pour 1 pop : 
data <- VCF1[, c("FORMAT",list_pop[[z]])]

# sliding window sur 1 pop 
test_100K <- sliding_spectre_onepop(data, 100000, 25000)
write_csv(test_100K, "/shared/projects/multiwhaling/multiwhaling/plot/scan_genom/SFS/test_100K.csv")

test_50K <- sliding_spectre_onepop(data, 50000, 10000)
write_csv(test_50K, "/shared/projects/multiwhaling/multiwhaling/plot/scan_genom/SFS/test_50K.csv")
