#### Script pour plot dans Admixture #### 
#
## input : pop_KX-combined-merged.txt ; issu de CLUMPP
## output : barplot d'appartenance des individus à différents cluster K
#


library(tidyverse)

# Charger données : 
admix <- read_table("prop_admix.txt", col_names =F)
popmap <- read_table("/shared/projects/multiwhaling/arivoltella/local_env/popmap.txt", col_names =F)
colnames(popmap) <- c("Ind", "pop")

# Pour ordre des individus sur le plot : 
order <- c("BERING", "KARAGINSKY", "CHILI", "PEROU", "MADAGASCAR", "N_ATL")


# Réarrangement du tibble : 
admix2 <- as_tibble(admix) |>
  dplyr::select(2:(ncol(admix)-1)) |>
  mutate(ind = popmap$Ind, 
         pop = factor(popmap$pop, levels = order)) |>
  arrange(pop) |>
  pivot_longer(cols = 1:(ncol(admix)-2), names_to = "K", values_to = "Ancestry_prop") |>
  mutate(K = gsub("X", "", K))

# Plot de structure : 
plot_admix <- admix2 |>
    mutate(ind = factor(ind, levels = unique(admix2$ind))) |>
  ggplot(aes(x = ind, y = Ancestry_prop, fill = K, color = K)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5, color = "black")) + 
  labs(x = "Individu", 
       y = "Ancestry proportions")

# Sauvegarde du graph : 
ggsave(paste("plot_admix_", (ncol(admix)-2), ".png", sep = ""), width = 7, height = 5)
saveRDS(admix2, paste("admix_", (ncol(admix)-2), ".RDS", sep = ""))


