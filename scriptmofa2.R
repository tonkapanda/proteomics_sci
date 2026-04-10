library(MOFA2)
library(tidyverse)


# section 1: load data

prot_raw <- read_tsv("human_proteomics.txt")
meta <- read_tsv("human_proteomics_meta.txt")

cat("proteins:", nrow(prot_raw), "| samples:", ncol(prot_raw) - 1, "\n")
cat("cohorts:", paste(unique(meta$Cohort), collapse = ", "), "\n")
cat("tissues:", paste(unique(meta$Tissue), collapse = ", "), "\n")
cat("timepoints:", paste(sort(unique(meta$Timepoint)), collapse = ", "), "\n")


# section 2: prepare matrices

prot_long <- prot_raw |>
  pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance") |>
  left_join(
    meta |> select(Sample, Tissue, Cohort, Timepoint, AIS_Base,
                   DeltaTMS, Conversion, Patient, Segment),
    by = "Sample"
  )

# exclude technical QC controls, keep patients and uninjured controls
discovery_samples <- meta |>
  filter(!str_detect(Sample, "PositiveControl|NegativeControl")) |>
  pull(Sample)

prot_discovery <- prot_long |>
  filter(Sample %in% discovery_samples)

# proteins x samples matrix for each biofluid
csf_mat <- prot_discovery |>
  filter(Tissue == "CSF") |>
  select(Gene, Sample, Abundance) |>
  pivot_wider(names_from = Sample, values_from = Abundance) |>
  column_to_rownames("Gene") |>
  as.matrix()

serum_mat <- prot_discovery |>
  filter(Tissue == "Serum") |>
  select(Gene, Sample, Abundance) |>
  pivot_wider(names_from = Sample, values_from = Abundance) |>
  column_to_rownames("Gene") |>
  as.matrix()

# keep proteins measured in at least 1/3 of samples, mirroring the paper
csf_mat <- csf_mat[rowMeans(!is.na(csf_mat)) >= 1/3, ]
serum_mat <- serum_mat[rowMeans(!is.na(serum_mat)) >= 1/3, ]

cat("\nafter missingness filter:\n")
cat("csf:", nrow(csf_mat), "proteins x", ncol(csf_mat), "samples\n")
cat("serum:", nrow(serum_mat), "proteins x", ncol(serum_mat), "samples\n")


# section 3: train models
# csf and serum have different sample sets so we run two separate models

# csf model
mofa_csf <- create_mofa(list(CSF = csf_mat))

data_opts_csf <- get_default_data_options(mofa_csf)
data_opts_csf$scale_views <- FALSE

model_opts_csf <- get_default_model_options(mofa_csf)
model_opts_csf$num_factors <- 10

train_opts_csf <- get_default_training_options(mofa_csf)
train_opts_csf$convergence_mode <- "slow"
train_opts_csf$seed <- 42

mofa_csf <- prepare_mofa(mofa_csf,
                         data_options = data_opts_csf,
                         model_options = model_opts_csf,
                         training_options = train_opts_csf)

mofa_csf_trained <- run_mofa(mofa_csf,
                             outfile = "MOFA2_CSF_model.hdf5",
                             use_basilisk = TRUE)

# serum model
mofa_serum <- create_mofa(list(Serum = serum_mat))

data_opts_serum <- get_default_data_options(mofa_serum)
data_opts_serum$scale_views <- FALSE

model_opts_serum <- get_default_model_options(mofa_serum)
model_opts_serum$num_factors <- 10

train_opts_serum <- get_default_training_options(mofa_serum)
train_opts_serum$convergence_mode <- "slow"
train_opts_serum$seed <- 42

mofa_serum <- prepare_mofa(mofa_serum,
                           data_options = data_opts_serum,
                           model_options = model_opts_serum,
                           training_options = train_opts_serum)

mofa_serum_trained <- run_mofa(mofa_serum,
                               outfile = "MOFA2_Serum_model.hdf5",
                               use_basilisk = TRUE)


# section 4: attach metadata to each trained model

meta_clean <- meta |>
  filter(!str_detect(Sample, "PositiveControl|NegativeControl")) |>
  select(Sample, Patient, Tissue, Timepoint, Cohort,
         AIS_Base, AIS_6mo, DeltaTMS, Conversion,
         Segment, Complete_6mo, AIS_Base_Numeric) |>
  mutate(
    Group = if_else(str_detect(Sample, "Control"), "Control", "SCI"),
    Timepoint_label = paste0(Timepoint, "h")
  ) |>
  rename(sample = Sample) |>
  as.data.frame()

samples_metadata(mofa_csf_trained) <- meta_clean |>
  filter(sample %in% unlist(samples_names(mofa_csf_trained)))

samples_metadata(mofa_serum_trained) <- meta_clean |>
  filter(sample %in% unlist(samples_names(mofa_serum_trained)))

# section 5: variance decomposition
# comparison point: paper found timepoint was the dominant source of variation

p_var_csf <- plot_variance_explained(mofa_csf_trained, plot_total = TRUE)
p_var_serum <- plot_variance_explained(mofa_serum_trained, plot_total = TRUE)

ggsave("MOFA2_variance_csf.pdf", p_var_csf, width = 8, height = 6)
ggsave("MOFA2_variance_serum.pdf", p_var_serum, width = 8, height = 6)


# section 6: factor scores coloured by clinical variables
# replicates the logic of figure 1d in the paper

p_csf_time <- plot_factors(mofa_csf_trained, factors = c(1, 2),
                           color_by = "Timepoint") +
  ggtitle("CSF: factor 1 vs 2 by timepoint")

p_csf_ais <- plot_factors(mofa_csf_trained, factors = c(1, 2),
                          color_by = "AIS_Base") +
  ggtitle("CSF: factor 1 vs 2 by AIS grade")

p_serum_time <- plot_factors(mofa_serum_trained, factors = c(1, 2),
                             color_by = "Timepoint") +
  ggtitle("Serum: factor 1 vs 2 by timepoint")

p_serum_ais <- plot_factors(mofa_serum_trained, factors = c(1, 2),
                            color_by = "AIS_Base") +
  ggtitle("Serum: factor 1 vs 2 by AIS grade")

ggsave("MOFA2_csf_factors_timepoint.pdf", p_csf_time, width = 8, height = 6)
ggsave("MOFA2_csf_factors_ais.pdf", p_csf_ais, width = 8, height = 6)
ggsave("MOFA2_serum_factors_timepoint.pdf", p_serum_time, width = 8, height = 6)
ggsave("MOFA2_serum_factors_ais.pdf", p_serum_ais, width = 8, height = 6)


# section 7: factor correlation with clinical covariates

correlate_factors_with_covariates(
  mofa_csf_trained,
  covariates = c("Timepoint", "AIS_Base_Numeric", "DeltaTMS",
                 "Complete_6mo", "Conversion"),
  plot = "log_pval"
)

correlate_factors_with_covariates(
  mofa_serum_trained,
  covariates = c("Timepoint", "AIS_Base_Numeric", "DeltaTMS",
                 "Complete_6mo", "Conversion"),
  plot = "log_pval"
)


# section 8: top protein weights per factor
# compare to paper's top csf biomarkers: GFAP, YWHAH, HSP90AA1, VIM, PRDX1
# serum biomarkers: LRG1, CRP, SERPINA5, AFM, APMAP

for (f in 1:5) {
  p <- plot_top_weights(mofa_csf_trained, view = "CSF",
                        factor = f, nfeatures = 15)
  ggsave(paste0("MOFA2_csf_weights_factor", f, ".pdf"), p, width = 8, height = 5)
  
  p <- plot_top_weights(mofa_serum_trained, view = "Serum",
                        factor = f, nfeatures = 15)
  ggsave(paste0("MOFA2_serum_weights_factor", f, ".pdf"), p, width = 8, height = 5)
}


# section 9: gfap weight profile
# gfap was the paper's most robust single biomarker in csf

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

ggsave("MOFA2_gfap_weights.pdf", p_gfap, width = 6, height = 4)




plot_top_weights(mofa_csf_trained, view = "CSF", factor = 1, nfeatures = 15)

plot_factors(mofa_csf_trained, factors = c(1,2), color_by = "Timepoint")

# section 10: save

saveRDS(mofa_csf_trained, "MOFA2_csf_trained.rds")
saveRDS(mofa_serum_trained, "MOFA2_serum_trained.rds")

cat("\ndone.\n")
sessionInfo()



