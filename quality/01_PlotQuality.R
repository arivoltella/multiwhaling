library(pcadapt)
library(vcfR)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(stringr)
library(LEA)
library(tidyverse)
library(hrbrthemes)
library(viridis)
#library(related)

## Read the VCF
VCF1 = read.vcfR("../Achille/data/sample_SNP.vcf.gz")

DP1 <- extract.gt(VCF1, element='DP', as.numeric = TRUE) 
DP_melt1 = melt(DP1)    # Pour avoir les noms d'individus 

geno1 <- as.data.frame(extract.gt(VCF1,element="GT", mask=F,as.numeric=F,return.alleles = F,
                 convertNA = F,extract = T))
# Génotype à chaque position pour chaque ind.


#### Populations #### -------------------------------------------------------------------
# Correspondance individus - population : 

names_ind <- colnames(DP1)
names_ind <- as_tibble(names_ind)
names_ind <- names_ind |>
  mutate(pop = str_extract(colnames(DP1), "^[A-Za-z]+|^[0-9]+")) |>
  mutate(pop = ifelse(pop == "A", "N_ATL",
               ifelse(pop == "Ber", "BERING", 
               ifelse(pop == "Chi", "CHILI", 
               ifelse(pop == "Kar", "KARAGINSKY", 
               ifelse(pop == "Pe", "PEROU", "MADAGASCAR"))))))

names_ind |>
  summarise(n = n(), 
            .by = pop) |>
  arrange(desc(n))


#### POSITION #### ----------------------------------------------------------------------
summary_position = data.frame(depth_pos = apply(DP1, 1, mean,na.rm=T))
position <- as.numeric(str_remove(rownames(summary_position), "1_"))

summary_position <- summary_position |>
  mutate(position = as.numeric(str_remove(rownames(summary_position), "1_")),  # Position du nucléotide 
         sum_NA = rowSums(geno1 == "./."),                         # Somme des NA / position
         prop_NA_ind = (sum_NA/dim(geno1)[2])*100,                 # Proportion des NA / individu
         het_pos = rowSums(geno1 == "0/1"),                        # Nb d'Hz / position
         prop_het = (het_pos/dim(geno1)[2])*100)                   # Proportion d'Hz / individu (%)
head(summary_position)




#### INDIVIDUS #### ---------------------------------------------------------------------
summary_individuals = data.frame(depth_ind = apply(DP1, 2, mean, na.rm = T))

summary_individuals <- summary_individuals |>
  mutate(ind = rownames(summary_individuals),                # Noms des individus 
         pop = names_ind$pop,                                # Population d'appartenance
         sum_NA = colSums(geno1=="./."),                     # Somme des NA / individu
         prop_NA_pos = ((sum_NA/dim(geno1)[1])*100),         # Proportion de NA / position
         het_ind = colSums(geno1=="0/1"),                    # Nb d'Hz / individu
         prop_het_pos = ((het_ind/dim(geno1)[1])*100))       # Proportion d'Hz / position
head(summary_individuals)




#### Plots of VCF before filtering #### -------------------------------------------------

#Depth distrib
g0 <- summary_position |>
  ggplot()+
  geom_histogram(aes(x = depth_pos),
		 binwidth = 4, fill = "steelblue", color = "black") +
  scale_x_continuous(limits = c(0, 200), breaks = seq(0, 100, by = 10)) +
  labs(x = "Distribition of average depth per site",
       y = "Count")


# Depth per site
g1 <- summary_position |>
  ggplot()+
  geom_point(aes(y = depth_pos, x = position),
             color="black", pch=20, size=1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y="mean depth", x = "PB")+
  ggtitle("Average depth sequencing of SNPs DP1")+
  scale_x_continuous(breaks = seq(0, 30e+06, by = 1000000)) +
  # vizualise only range of sequencing depth detph
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))


# Depth per Individual
g2 <- summary_individuals |>
  ggplot() + 
  geom_boxplot(aes(x = ind, y = depth_ind, color= pop)) +
  theme(axis.title.x = element_blank())+
  labs(y="SNP sequencing depth", x = "samples") +
  ggtitle("depth sequencing distribution of SNP filtered for qualtity and low depth") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # vizualise only range of sequencing depth detph
  scale_y_continuous(breaks = seq(0, 260, by = 20), limits = c(0, 260))


# NA per Position
g3 <- summary_position |>
  ggplot() +
  geom_point(aes(y = sum_NA, x = position), color="purple", size = 1)+
  labs(x="pos", y = "Prop './.' site") +
  ggtitle("Prop './.' site by position") +
  scale_x_continuous(breaks = seq(0, position[length(position)], 
                                  by = 1000000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# NA per Individual
g4 <- summary_individuals |>
  ggplot()+
  geom_point(aes(y = sum_NA, x = ind), color="purple", size=2) +
  labs(x="Ind", y = "Prop './.' site") +
  ggtitle("Prop './.' site by individuals") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))


# Het per Position
g5 <- summary_position |>
  ggplot() +
  geom_point(aes(y = het_pos , x= position), color = "purple", size = 0.3)+
  labs(x="pos", y = "Prop '0/1' site") +
  ggtitle("Prop '0/1' site by position") +
  scale_x_continuous(breaks = seq(0, position[length(position)], 
                                  by = 1000000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Het per Individual
g6 <- summary_individuals |>
  ggplot()+
  geom_point(aes(y = prop_het_pos, x = ind, color = pop), size = 1) +
  labs(x="Ind", y = "Prop '0/1' site") +
  ggtitle("Prop '0/1' site by individuals") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))



# Mean depth per position, per population : 
position = as.numeric(str_remove(rownames(DP1), "1_"))
DP2 <- as_tibble(DP1)
DP2 <- DP2 |>
  mutate(position = position) |>
  pivot_longer(cols = 1:88, 
               names_to = "individu", 
               values_to = "depth") |>
  mutate(pop = str_extract(individu, "^[A-Za-z]+|^[0-9]+")) |>
  mutate(pop = ifelse(pop == "A", "N_ATL",
                      ifelse(pop == "Ber", "BERING", 
                             ifelse(pop == "Chi", "CHILI", 
                                    ifelse(pop == "Kar", "KARAGINSKY", 
                                           ifelse(pop == "Pe", "PEROU", "MADAGASCAR"))))))
g7 <- DP2 |>
  group_by(position, pop) |>
  reframe(mean_DP = mean(depth)) |>
  ggplot(aes(x = position, y = mean_DP, color = pop)) + 
  geom_point(alpha = 0.3) +
  scale_y_continuous(limits = c(0,200), breaks = seq(0, 200, by = 10))


## Save plots
png(paste0("OldGenomeInfosVCF/", args[1],"_",args[2],"/DistribDepth_",args[1],"_",args[2],".png"))
print(g0)
dev.off()
png(paste0("OldGenomeInfosVCF/", args[1],"_",args[2],"/DepthPerSite_",args[1],"_",args[2],".png"))
print(g1)
dev.off()
png(paste0("OldGenomeInfosVCF/", args[1],"_",args[2],"/DepthPerInd_",args[1],"_",args[2],".png"))
print(g2)
dev.off()
png(paste0("OldGenomeInfosVCF/", args[1],"_",args[2],"/NAPerSite_",args[1],"_",args[2],".png"))
print(g3)
dev.off()
png(paste0("OldGenomeInfosVCF/", args[1],"_",args[2],"/NAPerInd_",args[1],"_",args[2],".png"))
print(g4)
dev.off()
png(paste0("OldGenomeInfosVCF/", args[1],"_",args[2],"/HetPerSite_",args[1],"_",args[2],".png"))
print(g5)
dev.off()
png(paste0("OldGenomeInfosVCF/", args[1],"_",args[2],"/HetPerInd_",args[1],"_",args[2],".png"))
print(g6)
dev.off()










  
  




















