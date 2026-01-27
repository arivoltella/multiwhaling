#################### 27/01/26 #####################
############### Compute the FIS ###################
############# Achille RIVOLTELLA ##################

# Libraries : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)


#### Read the VCF #  
VCF1 = readRDS("data/VCF_filtered.RDS")
list_pop <- readRDS("data/list_pop.RDS")


# ---------------------------------------------------------------------------------------------
######################### FIS ######################### ---------------------------------------
# ---------------------------------------------------------------------------------------------

#### Fonctions utiles pour la suite ---
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

genotypes <- extract.gt(VCF1, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
positions <- getPOS(VCF1)
snps_tot <- length(positions)
n_ind <- dim(genotypes)[2]

count_het <- apply(genotypes, 2, fun_geno_mod) |>
  as.numeric() |> matrix(ncol=n_ind, nrow=snps_tot, byrow=F)          # matrice numerique 0/1

count_allele <- apply(genotypes, 2, fun_geno_allele) |>
  as.numeric() |> matrix(ncol = n_ind, nrow = snps_tot, byrow=F)      # matrice numerique 0/1/2


#### Boucle permettant de calculer Hz attendue et observée : 
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

fis_by_pop <- all_FIS |>
  ggplot(aes(x = positions, y = FIS, color = pop))+
  geom_point()

png(paste0("plot/FIS/fis_by_pop.png"))
print(fis_by_pop)
dev.off()

mean_fis <- all_FIS |>
  filter(!is.na(FIS)) |>
  reframe(mean_FIS = round(mean(FIS), digits = 4),
          .by = pop) 
write_csv(mean_fis, "plot/FIS/mean_fis_by_pop.csv")
