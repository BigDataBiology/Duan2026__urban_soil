library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(cowplot)
library(scales)
library(ComplexHeatmap)
library(geosphere)
library(ggbreak)
library(ggforce)
library(ggpmisc)
library(Cairo)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else NULL
if (is.null(script_file)) {
  script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
}
analysis_candidates <- c(
  if (!is.null(script_file)) file.path(dirname(script_file), ".."),
  "analysis", ".", ".."
)
valid_analysis <- vapply(analysis_candidates, function(path) {
  dir.exists(file.path(path, "figure7")) &&
    dir.exists(file.path(path, "pre-calculated_data")) &&
    dir.exists(file.path(path, "..", "resource_generation"))
}, logical(1))
if (!any(valid_analysis)) stop("Could not locate the repository analysis directory")
setwd(normalizePath(analysis_candidates[which(valid_analysis)[1]]))
dir.create("figure7", showWarnings = FALSE)

gene_classes <- data.frame(rbind(
  c("aac2p",	"aminoglycoside acetyltransferase",	"aac(2')", "aminoglycoside"),
  c("aac3_1",	"aminoglycoside acetyltransferase",	"aac(3)", "aminoglycoside"), c("aac3_2",	"aminoglycoside acetyltransferase",	"aac(3)", "aminoglycoside"),
  c("aac6p_1",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"), c("aac6p_2",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"),
  c("aac6p_3",	"aminoglycoside acetyltransferase",	"aac(6')", "aminoglycoside"), c("aph2b",	"aminoglycoside phosphotransferase", "aph(2'')", "aminoglycoside"),
  c("aph3p",	"aminoglycoside phosphotransferase", "aph(3')", "aminoglycoside"), c("aph6p",	"aminoglycoside phosphotransferase", "aph(6)", "aminoglycoside"),
  c("class_a",	"beta-lactamase A",	"beta-lactamase A", "beta-lactamase"), c("class_b1_b2",	"beta-lactamase B",	"beta-lactamase B1-B2", "beta-lactamase"),
  c("class_b3",	"beta-lactamase B",	"beta-lactamase B3", "beta-lactamase"), c("class_c",	"beta-lactamase C",	"beta-lactamase C", "beta-lactamase"),
  c("class_d1",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"), c("class_d2",	"beta-lactamase D",	"beta-lactamase D", "beta-lactamase"),
  c("erm_1",	"macrolide erm", "erm", "macrolide"), c("erm_2",	"macrolide erm", "erm",	'macrolide'),
  c("mph",	"macrolide mph", "mph",	"macrolide"), c("qnr",	"qnr", "qnr", "quinolone"), c("tet_efflux",	"tetracycline efflux pump", "tet efflux", "tetracycline"),
  c("tet_enzyme",	"tetracycline inactivating enzymes", "tet enzyme", "tetracycline"), c("tet_rpg",	"tetracycline ribosomal protection genes", "tet RPG", "tetracycline")))

metadata <- readxl::read_excel("../resource_generation/12_ARGs_fARGene/supp1_sample_metadata.xlsx")
metadata <- metadata %>% mutate(Sample_id = gsub("ample", "", Sample_id)) 

centroids <- read.table("../resource_generation/12_ARGs_fARGene/contigs/orfs-amino-centroids.txt", quote="\"", comment.char="")
genes <- read.table("../resource_generation/12_ARGs_fARGene/contigs/closest_hits_orfs_amino.tsv", quote="\"", comment.char="")
genes <- genes %>% rename(orf = V1, centroid = V2) %>% select(orf, centroid)
genes <- genes %>% mutate(centroid = ifelse(orf %in% centroids$V1, orf, centroid))

genes$sample <- sapply(strsplit(genes$orf, split = "@@@"), function(x) x[1])
genes$hmm <- sapply(strsplit(genes$orf, split = "@@@"), function(x) x[2])
genes$c_hmm <- sapply(strsplit(genes$centroid, split = "@@@"), function(x) x[2])
genes$sequence <- sapply(strsplit(genes$orf, split = "@@@"), function(x) x[3])
genes$contig <- sapply(strsplit(genes$sequence, split = "_"), function(x) paste(x[1], x[2], x[3], sep = "_"))
genes <- genes %>% mutate(hmm = c_hmm)
genes$c_hmm <- NULL


######################################
######################################

## add metadata 


genes <- genes %>% mutate(description = gene_classes$X2[match(hmm, gene_classes$X1)],
                          class = gene_classes$X3[match(hmm, gene_classes$X1)],
                          hclass = gene_classes$X4[match(hmm, gene_classes$X1)],
                          City = metadata$City[match(sample, metadata$Sample_id)],
                          Location = metadata$Location[match(sample, metadata$Sample_id)],
                          Date = metadata$Date[match(sample, metadata$Sample_id)],
                          Longitude = metadata$Longitude[match(sample, metadata$Sample_id)],
                          Latitude = metadata$Latitude[match(sample, metadata$Sample_id)])


genes <- genes %>% 
  mutate(
    hclass = 
      factor(hclass, levels = c("aminoglycoside", "beta-lactamase", 
                                "macrolide", "tetracycline", "qnr")))



## remove sample s4 and s11

#genes <- genes %>% 
#  filter(!sample %in% c("s4", "s11"))


# font for plots 
general_size <- 10

# palette for plots  
pal_8 <- brewer.pal(8, "Dark2")
pal_8 <- pal_8[c(1,6,2:5,7,8)]
pal_div <- brewer.pal(8, "BrBG")
pal_seq <- brewer.pal(8, "YlOrBr")


## Calculate the pairwise sets
sets <- genes %>% 
  ungroup()  %>% 
  arrange(centroid) %>% 
  group_by(centroid) %>% 
  mutate(n_sample = n_distinct(sample)) %>% 
  mutate(single = (n_sample ==1))  %>%  
  group_by(sample) %>%
  summarise(centroids = list(centroid), .groups = "drop") %>% # put every query in a list
  mutate(Longitude = metadata$Longitude[match(sample, metadata$Sample_id)],
         Latitude = metadata$Latitude[match(sample, metadata$Sample_id)]) 

## calculate jaccard's index
JI <- expand_grid(sample_ref = sets$sample, sample_comp = sets$sample)  %>%
  left_join(sets, by = c("sample_ref" = "sample")) %>%
  rename(values1 = centroids) %>%
  left_join(sets, by = c("sample_comp" = "sample")) %>%
  rename(values2 = centroids) %>%
  mutate(jaccard = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length(union(.x, .y)))) %>% 
  mutate(n_overlap = map2_dbl(values1, values2, ~ length(intersect(.x, .y)))) %>% 
  mutate(n_ref = map2_dbl(values1, values2, ~ length(.x))) %>% 
  mutate(n_comp = map2_dbl(values1, values2, ~ length(.y))) %>% 
  group_by(sample_ref, sample_comp) %>% 
  mutate(distance = distHaversine(c(Longitude.x[1], Latitude.x), c(Longitude.y, Latitude.y))) %>%
  select(sample_ref, sample_comp, jaccard, distance, n_ref,  n_comp, n_overlap)  %>% 
  filter(sample_ref != sample_comp)


set.seed(2026)
JI <- JI[sample(1:dim(JI)[1], replace = F),]

## transform the jaccard index to matrix to use complex::heatmap
mat <- xtabs(jaccard ~ sample_ref + sample_comp, data = JI)
mat <- mat[,gtools::mixedsort(colnames(mat))]
mat <- mat[gtools::mixedsort(rownames(mat)),]
mat <- mat[metadata$Sample_id[metadata$Sample_id %in% genes$sample], 
           metadata$Sample_id[metadata$Sample_id %in% genes$sample]]

cities <- sort(unique(metadata$City))

# annotation for heatmap
col_ann = HeatmapAnnotation(
  City = metadata$City[match(colnames(mat), metadata$Sample_id) ],
  col = list(City = setNames(pal_8[1:2], cities)),
  annotation_name_side = NULL,
  show_annotation_name = FALSE)

row_ann = rowAnnotation(
  City     = metadata$City[match(rownames(mat), metadata$Sample_id)],
  col = list(City = setNames(pal_8[1:2], cities)),
  show_legend = FALSE,
  annotation_name_side = NULL,
  show_annotation_name = FALSE)



col_ann_nolegend = HeatmapAnnotation(
  City = metadata$City[match(colnames(mat), metadata$Sample_id) ],
  col = list(City = setNames(pal_8[1:2], cities)),
  annotation_name_side = NULL,
  show_annotation_name = FALSE,
  show_legend = FALSE)

row_ann_nolegend = rowAnnotation(
  City     = metadata$City[match(rownames(mat), metadata$Sample_id)],
  col = list(City = setNames(pal_8[1:2], cities)),
  show_legend = FALSE,
  annotation_name_side = NULL,
  show_annotation_name = FALSE)

# heatmap Jaccard's index
pheat11 <- ComplexHeatmap::Heatmap(
  mat, cluster_rows = T, cluster_columns = T,
  rect_gp = gpar(type = "none"), column_dend_side = "top",
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(as.numeric(x) >= 1 - as.numeric(y)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }},
  show_row_dend = FALSE,
  name = "J.I.", right_annotation = row_ann,
  top_annotation = col_ann, col = pal_seq, na_col = "white",
  show_row_names = FALSE, show_column_names = FALSE,
  heatmap_legend_param = list(direction = "horizontal"))


pheat12 <- ComplexHeatmap::Heatmap(
  mat, cluster_rows = T, cluster_columns = T,
  rect_gp = gpar(type = "none"), column_dend_side = "top",
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(as.numeric(x) >= 1 - as.numeric(y)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }},
  show_row_dend = FALSE,
  name = "J.I.", right_annotation = row_ann_nolegend,
  top_annotation = col_ann_nolegend, col = pal_seq, na_col = "white",
  show_row_names = FALSE, show_column_names = FALSE,
  heatmap_legend_param = list(direction = "horizontal"),
  show_heatmap_legend = FALSE)



##########################################################################################
##########################################################################################
##########################################################################################v


genes <- genes %>% 
  mutate(
    hclass = 
      factor(hclass, levels = c("aminoglycoside", "beta-lactamase", 
                                "macrolide", "tetracycline", "qnr")))


##########################################################################################
##########################################################################################
##########################################################################################


two_samples <- genes %>% select(centroid, City, Location, sample) %>% 
  group_by(centroid) %>% 
  mutate(N = n_distinct(sample)) %>% 
  filter(N > 1) %>% ungroup() %>% select(centroid) %>% distinct() %>% 
  pull() 

two_locations <- genes %>% select(centroid, City, Location, sample) %>% 
  group_by(centroid) %>% 
  mutate(N = n_distinct(Location)) %>% 
  filter(N > 1) %>% ungroup() %>% select(centroid) %>% distinct() %>% 
  pull() 

count_locations <- rbind(
  genes %>% 
    filter(centroid %in% two_locations) %>% 
    select(centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(City)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single city", "Multi city"),
           lev = "City"),
  
  genes %>% 
    filter(centroid %in% two_samples) %>% 
    select(centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(Location)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single location", "Multi location"),
           lev = "location"),
  
  
  genes %>% select(centroid, City, Location, sample) %>% 
    group_by(centroid) %>% 
    summarise(N = n_distinct(sample)) %>% 
    ungroup() %>% 
    group_by(N) %>% 
    summarise(n = n()) %>% 
    mutate(gr = ifelse(N == 1, "Single sample", "Multi sample"),
           lev = "sample"))


count_locations <- count_locations %>%
  ungroup() %>% 
  group_by(gr, lev) %>% 
  summarise(n = sum(n))


count_locations <- count_locations %>%
  mutate(lev = factor(lev, levels = c("mag", "sample", "location", "City"))) %>% 
  ungroup() %>% 
  group_by(lev) %>%
  mutate(p = n / sum(n))

count_locations <- count_locations %>% 
  mutate(gr2 = ifelse(grepl("Multi", gr), "2", "1"))


p_circle1 <- ggplot(count_locations, 
                    aes(x0 = 1 , y0 = 0, fill = gr2)) +
  geom_circle(aes(r = p), alpha = 1, show.legend = F) +
  scale_fill_manual(values = pal_8[3:8]) +
  geom_text(aes(x = 1 - p, y = p - 0.1, label = paste0(round((p) * 100), "%")), 
            size = (general_size-2) / ggplot2::.pt, vjust = 0, hjust = 0) +
  geom_text(aes(x = 1 - p, y = p - 0.1, label = paste0(n)), 
            size = (general_size-2) / ggplot2::.pt, vjust = 2.5, hjust = -0.5) +
  coord_equal() +
  theme_void() +
  labs(fill = "") + 
  guides(fill = guide_legend(ncol = 1)) +
  facet_wrap(lev ~ gr, ncol = 2) +
  theme(
    strip.text = element_text(size = general_size-2),
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    legend.text = element_text(size = general_size)) 

count_locations$gr2 <- factor(count_locations$gr2, levels = c(1,2))

p_circle1 <- ggplot(count_locations, aes(x = "", y = p, fill = gr2)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y", start = pi/1.5) +
  geom_text(aes(label = comma(n)),
            position = position_stack(vjust = 0.5),
            size = (general_size) / ggplot2::.pt) +
  facet_wrap(lev ~ ., ncol = 1) +
  theme_void() +
  scale_fill_manual(values = c(pal_8[4:3],pal_8[4:3],pal_8[4:3])) +
  theme(
    strip.text = element_text(size = general_size, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.text.y  = element_blank(),
    legend.text = element_text(size = general_size),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank())


# Shared input for Fig. 7a; the notebook owns its plot and whisker definition.
table_plot <- genes %>%
  group_by(City, Location, sample) %>%
  summarise(n = n_distinct(centroid), .groups = "drop")

write.csv(table_plot, "pre-calculated_data/number_ARGs_per_sample.csv", row.names = FALSE)

# Fig. 7b: sharing of ARG centroids across samples, locations and cities.
ggsave("figure7/7_b.svg", p_circle1 +
         theme(plot.margin = margin(0, 0, 0, 0), 
               panel.spacing = unit(5, "mm")) +
         ggtitle("B"), 
       width = 90, height = 100, unit = "mm")


# Fig. 7c: pairwise Jaccard similarity.
svg("figure7/7_c.svg", width = 70/25.4, height = 100/25.4)
draw(pheat11, heatmap_legend_side = "bottom", merge_legend = TRUE)
grid.text("C", x = 0.5, y = unit(1, "npc") + unit(2, "mm"), gp = gpar(fontsize=general_size + 2))
dev.off()




# Auxiliary analyses below do not draw Fig. 7a.
per_sample_centroids <- lapply(split(genes$centroid, genes$sample), unique)
all_samples <- names(per_sample_centroids)
n <- length(all_samples)
max_j <- min(58, n)

grid <- expand.grid(j = seq_len(max_j), k = seq_len(50))

out <- Map(function(j, k) {
  chosen <- all_samples[sample.int(n, size = j, replace = FALSE)]
  c(k = k, j = j, tot = length(unique(unlist(per_sample_centroids[chosen]))))
}, grid$j, grid$k)

results_genes <- as.data.frame(do.call(rbind, out))

summary_genes <- results_genes %>%
  group_by(j) %>%
  summarise(mean_tot = mean(tot), sd_tot = sd(tot), .groups = "drop")

ggplot(summary_genes, aes(x = j, y = mean_tot)) +
  geom_ribbon(aes(ymin = mean_tot - sd_tot, ymax = mean_tot + sd_tot), alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(x = "Number of samples (j)", y = "Unique centroids") +
  theme_minimal()

write.csv(results_genes %>% mutate(dataset = "UrbanSoil"), file = "../resource_generation/12_ARGs_fARGene/incremental_curve_soil.csv", row.names = F)

########




library(uwot)
library(Matrix)
library(vegan)

build_composition_matrix <- function(gene_tbl, weighted = FALSE) {
  if (weighted) {
    tbl <- gene_tbl %>%
      count(sample, centroid, name = "n_hits") %>%
      pivot_wider(names_from = centroid, values_from = n_hits, values_fill = 0)
  } else {
    tbl <- gene_tbl %>%
      distinct(sample, centroid) %>%
      mutate(present = 1) %>%
      pivot_wider(names_from = centroid, values_from = present, values_fill = 0)
  }
  mat <- as.matrix(tbl[, -1])
  rownames(mat) <- tbl$sample
  mat
}

pa_mat    <- build_composition_matrix(genes, weighted = FALSE)

set.seed(2026)
jacc_dist <- vegdist(pa_mat, method = "jaccard", binary = TRUE)
pcoa <- cmdscale(jacc_dist, k = 2, eig = TRUE, add = TRUE)
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCo1", "PCo2")
pcoa_df$sample   <- rownames(pa_mat)
pcoa_df$City     <- metadata$City[match(pcoa_df$sample, metadata$Sample_id)]
pcoa_df$Location <- metadata$Location[match(pcoa_df$sample, metadata$Sample_id)]

eig  <- pcoa$eig
pct1 <- round(100 * eig[1] / sum(eig[eig > 0]), 1)
pct2 <- round(100 * eig[2] / sum(eig[eig > 0]), 1)


#locations  <- sort(unique(metadata$Location[match(sample_order, metadata$Sample_id)]))

pcoa_plot <- ggplot(pcoa_df, aes(PCo1, PCo2, color = City, shape = Location)) +
  geom_point(size = 3) +
  scale_color_manual(values = setNames(pal_8[1:2], cities)) +
  #scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8)[seq_along(locations)]) +
  labs(x = paste0("PCo1 (", pct1, "%)"), y = paste0("PCo2 (", pct2, "%)"),
       title = "PCoA (Jaccard)") +
  theme_minimal()

ggsave("figure7/pcoa_args.svg", pcoa_plot +
         theme(plot.margin = margin(0, 0, 0, 0), 
               panel.spacing = unit(5, "mm")),
       width = 140, height = 100, unit = "mm")


write.csv( pcoa_df, "figure7/PCoA_ARGs.csv")
