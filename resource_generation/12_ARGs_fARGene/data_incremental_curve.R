

unigenes <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/unigenes_per_tool.rds")
args_abundances <- read.delim("data/abundances/args_abundances.tsv")
metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
metadata <- metadata %>% mutate(sample = sample_id)
args_abundances <- args_abundances %>% left_join(metadata[,c("sample","insertsHQ", "insertsRaw")], by = "sample")

genes_right_habitat <- c("human gut", "wastewater", "marine", "freshwater", "soil" )

fg <- unigenes$query[unigenes$tool %in% "fARGene"]

args_abundances_sb <- args_abundances %>% 
  filter(X %in% fg) %>%
  mutate(habitat = metadata$habitat[match(sample, metadata$sample)]) %>% 
  filter(habitat %in% genes_right_habitat, sample  %in% metadata$sample_id) 


set.seed(2025)

rarefaction_fast <- function(X, raw, inserts, depth = 5e6, seed = 2025){
  # set.seed(seed)  # optional
  reads_count <- ceiling(raw)
  arg_reads <- sum(reads_count)
  not_arg_reads <- inserts[1] - arg_reads
  # adjust depth
  depth <- min(depth, inserts[1])
  # probabilities
  probs <- c(reads_count, not_arg_reads) / inserts[1]
  names(probs) <- c(X, "not an ARG")
  # multinomial sampling
  sampled <- rmultinom(1, size = depth, prob = probs)
  return(sampled[,1])
}


rarefied_counts <- args_abundances_sb %>%
  filter(raw > 0) %>%
  ungroup() %>%
  group_by(sample) %>%
  summarise(
    rarefied = list(
      rarefaction_fast(X, raw, insertsHQ, 5e6)
    ),
    .groups = "drop"
  ) %>%
  unnest_longer(rarefied, indices_include = TRUE) %>%
  rename(X = rarefied_id, rarified_count = rarefied)


# add them to the abundance object   
args_abundances_sb <- args_abundances_sb %>% left_join(rarefied_counts, by = c("sample", "X"))
args_abundances_sb <- args_abundances_sb %>% mutate(rarified_count = ifelse(is.na(rarified_count), 0, rarified_count))



# Precompute once: for each sample, its set of ARGs (total, and rarified-only)
per_sample <- lapply(split(args_abundances_sb, args_abundances_sb$sample), function(d) {
  list(all = unique(d$X), rar = unique(d$X[d$rarified_count > 0]), habitat = d$habitat[1])
})

samples_by_habitat <- split(names(per_sample), sapply(per_sample, `[[`, "habitat"))
rarefy_habitat <- function(samples_h, habitat_name, per_sample, n_reps = 50, max_size = 100) {
  n <- length(samples_h)
  max_j <- min(max_size, n)
  grid <- expand.grid(j = seq_len(max_j), k = seq_len(n_reps))
  
  out <- Map(function(j, k) {
    chosen <- samples_h[sample.int(n, size = j, replace = FALSE)]
    c(k = k, j = j,
      tot = length(unique(unlist(lapply(per_sample[chosen], `[[`, "all")))),
      rar = length(unique(unlist(lapply(per_sample[chosen], `[[`, "rar")))))
  }, grid$j, grid$k)
  
  data.frame(do.call(rbind, out), habitat = habitat_name)
}

results <- do.call(rbind, Map(rarefy_habitat, samples_by_habitat, names(samples_by_habitat),
                              MoreArgs = list(per_sample = per_sample)))

summary_df <- results %>%
  group_by(habitat, j) %>%
  summarise(mean_tot = mean(tot), sd_tot = sd(tot), .groups = "drop")

ggplot(summary_df, aes(x = j, y = mean_tot, color = habitat, fill = habitat)) +
  geom_ribbon(aes(ymin = mean_tot - sd_tot, ymax = mean_tot + sd_tot), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Number of samples (j)", y = "Unique predicted ARGs", color = "Habitat", fill = "Habitat") +
  theme_minimal()


write.csv(results %>% mutate(dataset = "GMGC"), file = "~/Documents/GitHub/urban_soil/resource_generation/12_ARGs_fARGene/incremental_curve_gmgc.csv", row.names = F)

