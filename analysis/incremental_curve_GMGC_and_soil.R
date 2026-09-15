
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(grid)
library(cowplot)
library(ggbreak)
library(ggforce)
library(ggpmisc)

gmgc <- read.csv("~/Documents/GitHub/urban_soil/resource_generation/12_ARGs_fARGene/incremental_curve_gmgc.csv")
urbansoil <- read.csv("~/Documents/GitHub/urban_soil/resource_generation/12_ARGs_fARGene/incremental_curve_soil.csv")
urbansoil_sr <- read.csv("~/Documents/GitHub/urban_soil/resource_generation/12_ARGs_fARGene/incremental_curve_soil_sr.csv")

summary_df <- bind_rows(gmgc , urbansoil %>% mutate(habitat = "Urban Soil - LR"), urbansoil_sr %>% mutate(habitat = "Urban Soil - SR")) %>%
  group_by(habitat, j) %>%
  summarise(mean_tot = mean(tot), sd_tot = sd(tot), .groups = "drop")

g <- ggplot(summary_df, aes(x = j, y = mean_tot, color = habitat, fill = habitat)) +
  geom_ribbon(aes(ymin = mean_tot - sd_tot, ymax = mean_tot + sd_tot), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Number of samples", y = "Unique predicted ARGs", color = "Habitat", fill = "Habitat") +
  theme_minimal()


svg("supp9/incremental_curve.svg", width = 10, height = 8)
g
dev.off()


pdf("supp9/incremental_curve.pdf", width = 10, height = 8)
g
dev.off()


data.frame(summary_df[summary_df$j<60 & (!summary_df$habitat  %in% c("freshwater","marine")),])
