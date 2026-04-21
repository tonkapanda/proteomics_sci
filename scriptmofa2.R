# 1. setup
library(MOFA2)
library(tidyverse)

dir.create("results", showWarnings = FALSE)

# 2. load data
prot_raw <- read_tsv("human_proteomics.txt")
meta     <- read_tsv("human_proteomics_meta.txt")

cat("proteins:", nrow(prot_raw), "| samples:", ncol(prot_raw) - 1, "\n")
cat("cohorts:",  paste(unique(meta$Cohort),    collapse = ", "), "\n")
cat("tissues:",  paste(unique(meta$Tissue),    collapse = ", "), "\n")
cat("timepoints:", paste(sort(unique(meta$Timepoint)), collapse = ", "), "\n")

# 3. long-format join and qc filter
prot_long <- prot_raw |>
  pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance") |>
  left_join(
    meta |> select(Sample, Tissue, Cohort, Timepoint, AIS_Base,
                   DeltaTMS, Conversion, Patient, Segment),
    by = "Sample"
  )

# drop technical qc controls
discovery_samples <- meta |>
  filter(!str_detect(Sample, "PositiveControl|NegativeControl")) |>
  pull(Sample)

prot_discovery <- prot_long |>
  filter(Sample %in% discovery_samples)

# 4. build per-biofluid matrices
# helper: pivot to proteins x samples matrix
make_mat <- function(df, tissue) {
  df |>
    filter(Tissue == tissue) |>
    select(Gene, Sample, Abundance) |>
    pivot_wider(names_from = Sample, values_from = Abundance) |>
    column_to_rownames("Gene") |>
    as.matrix()
}

csf_mat   <- make_mat(prot_discovery, "CSF")
serum_mat <- make_mat(prot_discovery, "Serum")

# require protein measured in >= 1/3 of samples (mirrors paper)
csf_mat   <- csf_mat[rowMeans(!is.na(csf_mat))   >= 1/3, ]
serum_mat <- serum_mat[rowMeans(!is.na(serum_mat)) >= 1/3, ]

cat("\nafter missingness filter:\n")
cat("csf:",   nrow(csf_mat),   "proteins x", ncol(csf_mat),   "samples\n")
cat("serum:", nrow(serum_mat), "proteins x", ncol(serum_mat), "samples\n")

# 5. train mofa2 models
# helper: configure and run a single-view mofa model
train_mofa <- function(mat, view_name, outfile, k = 10, seed = 42) {
  obj <- create_mofa(setNames(list(mat), view_name))

  d_opts <- get_default_data_options(obj)
  d_opts$scale_views <- FALSE

  m_opts <- get_default_model_options(obj)
  m_opts$num_factors <- k

  t_opts <- get_default_training_options(obj)
  t_opts$convergence_mode <- "slow"
  t_opts$seed <- seed

  obj <- prepare_mofa(obj,
                      data_options     = d_opts,
                      model_options    = m_opts,
                      training_options = t_opts)

  run_mofa(obj, outfile = outfile, use_basilisk = TRUE)
}

# csf and serum have disjoint sample sets, so separate models
mofa_csf_trained   <- train_mofa(csf_mat,   "CSF",   "results/MOFA2_CSF_model.hdf5")
mofa_serum_trained <- train_mofa(serum_mat, "Serum", "results/MOFA2_Serum_model.hdf5")

# 6. sensitivity test on number of factors (k = 3-10)
k_range <- 3:10
sens_results <- list()

for (k in k_range) {
  cat("training csf  k =", k, "...\n")
  m_csf <- train_mofa(csf_mat, "CSF",
                       paste0("results/MOFA2_CSF_k", k, ".hdf5"), k = k)
  r2_csf <- calculate_variance_explained(m_csf)$r2_total[[1]]

  cat("training serum k =", k, "...\n")
  m_ser <- train_mofa(serum_mat, "Serum",
                       paste0("results/MOFA2_Serum_k", k, ".hdf5"), k = k)
  r2_ser <- calculate_variance_explained(m_ser)$r2_total[[1]]

  sens_results[[length(sens_results) + 1]] <- tibble(
    k = k, tissue = "CSF", total_r2 = r2_csf
  )
  sens_results[[length(sens_results) + 1]] <- tibble(
    k = k, tissue = "Serum", total_r2 = r2_ser
  )
}

sens_df <- bind_rows(sens_results)

# total variance explained vs number of factors
p_sens <- ggplot(sens_df, aes(x = k, y = total_r2, colour = tissue)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = k_range) +
  labs(title = "MOFA2 sensitivity: total variance explained by number of factors",
       x = "number of factors (k)", y = "total variance explained (R²)",
       colour = "biofluid") +
  theme_bw()

ggsave("results/MOFA2_factor_sensitivity.pdf", p_sens, width = 8, height = 5)
cat("\nsensitivity plot saved\n")

# 7. attach metadata to trained models
meta_clean <- meta |>
  filter(!str_detect(Sample, "PositiveControl|NegativeControl")) |>
  select(Sample, Patient, Tissue, Timepoint, Cohort,
         AIS_Base, AIS_6mo, DeltaTMS, Conversion,
         Segment, Complete_6mo, AIS_Base_Numeric) |>
  mutate(
    Group           = if_else(str_detect(Sample, "Control"), "Control", "SCI"),
    Timepoint_label = paste0(Timepoint, "h")
  ) |>
  rename(sample = Sample) |>
  as.data.frame()

samples_metadata(mofa_csf_trained) <- meta_clean |>
  filter(sample %in% unlist(samples_names(mofa_csf_trained)))

samples_metadata(mofa_serum_trained) <- meta_clean |>
  filter(sample %in% unlist(samples_names(mofa_serum_trained)))

# 8. variance decomposition
# paper found timepoint was the dominant source of variation
p_var_csf   <- plot_variance_explained(mofa_csf_trained,   plot_total = TRUE)
p_var_serum <- plot_variance_explained(mofa_serum_trained, plot_total = TRUE)

ggsave("results/MOFA2_variance_csf.pdf",   p_var_csf,   width = 8, height = 6)
ggsave("results/MOFA2_variance_serum.pdf", p_var_serum, width = 8, height = 6)

# 9. factor scatter plots (cf. fig 1d)
p_csf_time  <- plot_factors(mofa_csf_trained,   factors = c(1, 2), color_by = "Timepoint") +
  ggtitle("CSF: factor 1 vs 2 by timepoint")
p_csf_ais   <- plot_factors(mofa_csf_trained,   factors = c(1, 2), color_by = "AIS_Base") +
  ggtitle("CSF: factor 1 vs 2 by AIS grade")
p_serum_time <- plot_factors(mofa_serum_trained, factors = c(1, 2), color_by = "Timepoint") +
  ggtitle("Serum: factor 1 vs 2 by timepoint")
p_serum_ais  <- plot_factors(mofa_serum_trained, factors = c(1, 2), color_by = "AIS_Base") +
  ggtitle("Serum: factor 1 vs 2 by AIS grade")

ggsave("results/MOFA2_csf_factors_timepoint.pdf",   p_csf_time,  width = 8, height = 6)
ggsave("results/MOFA2_csf_factors_ais.pdf",         p_csf_ais,   width = 8, height = 6)
ggsave("results/MOFA2_serum_factors_timepoint.pdf", p_serum_time, width = 8, height = 6)
ggsave("results/MOFA2_serum_factors_ais.pdf",       p_serum_ais,  width = 8, height = 6)

# 10. factor-covariate correlations
cov_vars <- c("Timepoint", "AIS_Base_Numeric", "DeltaTMS",
              "Complete_6mo", "Conversion")

correlate_factors_with_covariates(mofa_csf_trained,   covariates = cov_vars, plot = "log_pval")
correlate_factors_with_covariates(mofa_serum_trained, covariates = cov_vars, plot = "log_pval")

# 11. top protein weights per factor
# paper's top csf: GFAP, YWHAH, HSP90AA1, VIM, PRDX1
# paper's top serum: LRG1, CRP, SERPINA5, AFM, APMAP
for (f in 1:5) {
  p <- plot_top_weights(mofa_csf_trained,   view = "CSF",   factor = f, nfeatures = 15)
  ggsave(paste0("results/MOFA2_csf_weights_factor", f, ".pdf"),   p, width = 8, height = 5)

  p <- plot_top_weights(mofa_serum_trained, view = "Serum", factor = f, nfeatures = 15)
  ggsave(paste0("results/MOFA2_serum_weights_factor", f, ".pdf"), p, width = 8, height = 5)
}

# 12. gfap deep dive
# gfap was the paper's most robust single csf biomarker
csf_weights <- get_weights(mofa_csf_trained, views = "CSF", as.data.frame = TRUE)

gfap_weights <- csf_weights |>
  filter(str_detect(feature, "GFAP")) |>
  arrange(desc(abs(value)))

cat("\nGFAP weights across factors (CSF):\n")
print(gfap_weights)

p_gfap <- gfap_weights |>
  ggplot(aes(x = factor, y = value)) +
  geom_col(fill = "steelblue") +
  labs(title = "GFAP weight across MOFA2 factors (CSF)",
       x = "factor", y = "weight") +
  theme_bw()

ggsave("results/MOFA2_gfap_weights.pdf", p_gfap, width = 6, height = 4)

# 13. save and session info
saveRDS(mofa_csf_trained,   "results/MOFA2_csf_trained.rds")
saveRDS(mofa_serum_trained, "results/MOFA2_serum_trained.rds")

cat("\ndone.\n")
sessionInfo()