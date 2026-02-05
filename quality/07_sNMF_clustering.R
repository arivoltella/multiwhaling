################### 27/01/26 ######################
########## sNMF clustering analysis ###############
############# Achille RIVOLTELLA ##################

# Library : 
library(vcfR)
library(rlist)
library(qqman)
library(ggplot2)
library(tidyverse)
library(pegas)
library(MASS)
library(patchwork)
library(LEA)


#### Read the data #  
VCF1 = readRDS("/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered.RDS")
list_pop <- readRDS("/shared/projects/multiwhaling/multiwhaling/data/list_pop.RDS")

names_ind <- enframe(list_pop, name = "Population", value = "Individu") |>
  unnest(Individu)

# ------------------------------------------------------------------------------------------
################## sNMF clustering ############## ------------------------------------------
# ------------------------------------------------------------------------------------------

whale_snmf <- snmf(input.file= "/shared/projects/multiwhaling/multiwhaling/data/VCF_filtered.geno", 
                   K= 1:6, repetitions = 3, project="new", entropy=T)
#show(whale_snmf)
summary(whale_snmf)

pdf(file = "/shared/projects/multiwhaling/multiwhaling/plot/Clustering/entropy_K.pdf")
plot(whale_snmf, col = "blue", pch = 19, cex = 1.2) ### Identifier le K avec l'entropie la plus faible
dev.off()

# Initialisation des objets : 
list_plot <- c()
order <- c("BERING", "KARAGINSKY", "CHILI", "PEROU", "MADAGASCAR", "N_ATL")

for (i in 3:6) {
  best <- which.min(cross.entropy(whale_snmf, K = i))
  adm_coeff <- Q(whale_snmf, K = i, run = best)
  
  adm_coeff <- as_tibble(adm_coeff) |>
    mutate(ind = names_ind$Individu, 
           pop = factor(names_ind$Population, levels = order)) |>
    arrange(pop) |>
    pivot_longer(cols = 1:i, names_to = "K", values_to = "Ancestry_prop") |>
    mutate(K = gsub("V", "", K))
  
  plot_structure <- adm_coeff |>
    mutate(ind = factor(ind, levels = unique(adm_coeff$ind))) |>
    ggplot(aes(x = ind, y = Ancestry_prop, fill = K, color = K)) + 
    geom_col() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 3, color = "black")) + 
    labs(x = "Individu", 
         y = "Ancestry proportions")
  
  list_plot <- c(list_plot, list(plot_structure))
}

structure <- wrap_plots(list_plot) + plot_layout(axes = "collect")
ggsave("/shared/projects/multiwhaling/multiwhaling/plot/Clustering/structure.png", height = 5, width = 9)



# best <- which.min(cross.entropy(whale_snmf, K = 3))
# adm_coeff <- Q(whale_snmf, K = 3, run = best)
# pdf(file = "plot/Clustering/structure_plot3.pdf")
# barplot(t(adm_coeff), names = names_ind$Individu, las = 2, cex.names = 0.6, border = NA, space = 0,
#         col = my.colors,
#         xlab = "Individuals",
#         ylab = "Ancestry proportions",
#         main = "Ancestry matrix")
# dev.off()