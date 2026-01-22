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

args <- commandArgs(trailingOnly = TRUE)

## Read the VCF

VCF1=read.vcfR(paste0("VCFNewNamesOldGenome/",args[2],"_NewGenome_",args[1],"_GATK_TAG_Flowqual_Noindels_Norepeat_SNP.vcf.gz"))

DP1<- extract.gt(VCF1, element='DP', as.numeric = TRUE) 
DP_melt1=melt(DP1)
geno1<-extract.gt(VCF1,element="GT",mask=F,as.numeric=F,return.alleles = F,
                 convertNA = F,extract = T)
geno_table1=as.data.frame(geno1)

#Coverage Information

moyenne_cover_par_indiv1 = apply(DP1,2,mean,na.rm=T)
moyenne_cover_par_Ind_table1=as.data.frame(moyenne_cover_par_indiv1)

#Mean Per Line

moyenne_cover_par_pos1= apply(DP1,1,mean,na.rm=T)
moyenne_cover_par_pos_table1=as.data.frame(moyenne_cover_par_pos1)

#Position vector

position_vector1=as.numeric(str_remove(rownames(DP1),paste0(args[1],"_")))
moyenne_cover_par_pos_table_pos1=cbind(moyenne_cover_par_pos_table1, position_vector1)

#Row count of Genotype NA

sum_pos_geno_table1 = as.data.frame(rowSums(geno1=="./.")) # put "./.", "0/0", "0/1" or "1/1" depending on the genotype you want to count

# By individuals

sum_ind_geno_table1 = as.data.frame(colSums(geno1=="./."))

# By position

Prop_Pos_geno1=((sum_pos_geno_table1/dim(geno_table1)[2])*100)
Prop_Pos_geno_table1=as.data.frame(Prop_Pos_geno1)
position_vector1=as.numeric(str_remove(rownames(Prop_Pos_geno_table1),paste0(args[1],"_")))
Prop_Pos_geno_table_pos1=cbind(Prop_Pos_geno_table1, position_vector1)

# By individuals

Prop_Ind_geno1=((sum_ind_geno_table1/dim(geno_table1)[1])*100)
Prop_Ind_geno_table1=as.data.frame(Prop_Ind_geno1)

#Row count of Het

sum_pos_geno_table1_Het = as.data.frame(rowSums(geno1=="0/1")) # put "./.", "0/0", "0/1" or "1/1" depending on the genotype you want to count

# By individuals

sum_ind_geno_table1_Het = as.data.frame(colSums(geno1=="0/1"))

# By position

Prop_Pos_geno1_Het=((sum_pos_geno_table1_Het/dim(geno_table1)[2])*100)
Prop_Pos_geno_table1_Het=as.data.frame(Prop_Pos_geno1_Het)
Prop_Pos_geno_table_pos1_Het=cbind(Prop_Pos_geno_table1_Het, position_vector1)

# By individuals
Prop_Ind_geno1_Het=((sum_ind_geno_table1_Het/dim(geno_table1)[1])*100)
Prop_Ind_geno_table1_Het=as.data.frame(Prop_Ind_geno1_Het)

## Plots of VCF before filtering

#Depth distrib

g0=ggplot()+
  geom_histogram(data=moyenne_cover_par_pos_table_pos1,
		 aes(x=moyenne_cover_par_pos_table_pos1$moyenne_cover_par_pos1),
		 binwidth = 4, fill = "steelblue", color = "black") +
  scale_x_continuous(limits = c(0, 200), breaks = seq(0, 100, by = 10)) +
  labs(x =paste0("Distribition of average depth per site (",args[2]," Chrom ",args[1],")"),
       y = "Count")
# Depth per site

g1=ggplot()+
  geom_point(data=moyenne_cover_par_pos_table_pos1, 
             aes(y=moyenne_cover_par_pos_table_pos1$moyenne_cover_par_pos1,
                 x=moyenne_cover_par_pos_table_pos1$position_vector1),
             color="black", pch=20, size=1) +
  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y="mean depth", x = "PB")+
  ggtitle("Average depth sequencing of SNPs DP1")+
  scale_x_continuous(breaks = seq(0, 30e+06, by = 1000000)) +
  # vizualise only range of sequencing depth detph
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))

# Depth per Individual

g2=ggplot()+
geom_boxplot(data=DP_melt1, aes(x=DP_melt1$Var2, 
                               y=DP_melt1$value,
                                 colour="purple"),
               color="black") +
  theme(axis.title.x=element_blank())+
  labs(y="SNP sequencing depth", x = "samples")+
  ggtitle("depth sequencing distribution of SNP filtered for qualtity and low depth") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  # vizualise only range of sequencing depth detph
  scale_y_continuous(breaks = seq(0, 260, by = 20), limits = c(0, 260))

# NA per Position

g3=ggplot()+
  geom_point(data=Prop_Pos_geno_table_pos1, 
             aes(y=Prop_Pos_geno_table_pos1$`rowSums(geno1 == "./.")` , 
                 x=Prop_Pos_geno_table_pos1$position_vector1), 
             color="purple",
             size=1)+
  labs(x="pos", y = "Prop './.' site") +
  ggtitle("Prop './.' site by position") +
  scale_x_continuous(breaks = seq(0, position_vector1[length(position_vector1)], 
                                  by = 1000000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# NA per Individual

g4=ggplot()+
  geom_point(data=Prop_Ind_geno_table1, 
             aes(y=Prop_Ind_geno_table1$`colSums(geno1 == "./.")`, 
                 x=row.names(Prop_Ind_geno_table1)),
             color="purple",
             size=2) +
  labs(x="Ind", y = "Prop './.' site") +
  ggtitle("Prop './.' site by individuals") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))

# Het per Position

g5=ggplot()+
  geom_point(data=Prop_Pos_geno_table_pos1_Het, 
             aes(y=Prop_Pos_geno_table_pos1_Het$`rowSums(geno1 == "0/1")` , 
                 x=Prop_Pos_geno_table_pos1_Het$position_vector1), 
             color="purple",
             size=1)+
  labs(x="pos", y = "Prop '0/1' site") +
  ggtitle("Prop '0/1' site by position") +
  scale_x_continuous(breaks = seq(0, position_vector1[length(position_vector1)], 
                                  by = 1000000)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Het per Individual

g6=ggplot()+
  geom_point(data=Prop_Ind_geno_table1_Het, 
             aes(y=Prop_Ind_geno_table1_Het$`colSums(geno1 == "0/1")`, 
                 x=row.names(Prop_Ind_geno_table1_Het)),
             color="purple",
             size=2) +
  labs(x="Ind", y = "Prop '0/1' site") +
  ggtitle("Prop '0/1' site by individuals") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100))

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
