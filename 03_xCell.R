# Please run 02_tpm.R, folllowed by xCell before running this script.

# Load libraries
library(tidyverse)
library(sva)
library(ggplot2)
library(gridExtra) # Useful for side-by-side plotting

# Create output directories
dir.create("./results/xCell_results/Corrected_scores/", showWarnings = FALSE, recursive = TRUE)

# Define Color Scheme for PCA plots
my_colors <- c(
  "Control" = "#79AC78",
  "HP" = "#D25353",
  "IPF" = "#5D688A",
  "GSE184316" = "#EF9C66",
  "GSE150910" = "#9ECAD6"
)

# Read the xCell cell type enrichment scores
# "es" means Enrichment scores
gse184316_hp_ctrl_es <- read.table("./results/xCell_results/GSE184316_HP_Ctrl/xCell_TPM_GSE184316_HP_Ctrls_xCell_0656021926.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

gse184316_hp_ipf_es <- read.table("./results/xCell_results/GSE184316_HP_IPF/xCell_TPM_GSE184316_HP_IPF_xCell_0659021926.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

gse184316_ipf_ctrl_es <- read.table("./results/xCell_results/GSE184316_IPF_Ctrl/xCell_TPM_GSE184316_IPF_Ctrls_xCell_0703021926.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

gse150910_hp_ctrl_es <- read.table("./results/xCell_results/GSE150910_HP_Ctrl/xCell_TPM_GSE150910_HP_Ctrls_xCell_0708021926.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

gse150910_hp_ipf_es <- read.table("./results/xCell_results/GSE150910_HP_IPF/xCell_TPM_GSE150910_HP_IPF_xCell_0712021926.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

gse150910_ipf_ctrl_es <- read.table("./results/xCell_results/GSE150910_IPF_Ctrl/xCell_TPM_GSE150910_IPF_Ctrls_xCell_0718021926.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)


# Now we will run PCA by merging each condition together
# eg: We will merge GSE184316_HP_Ctrl with GSE150910_HP_Ctrl

## Condition 1: cHP vs Control #############################
# First we will merge metadata
gse184316_metadata_chp_ctrl <- read.csv("./data/GSE184316/Sample_information_HP_Control_GSE184316.csv",
  header = TRUE,
  row.names = 1
)

gse150910_metadata_chp_ctrl <- read.csv("./data/GSE150910/Sample_information_HP_Control_GSE150910.csv",
  header = TRUE,
  row.names = 1
)

merged_metadata_chp_ctrl <- rbind(
  gse184316_metadata_chp_ctrl,
  gse150910_metadata_chp_ctrl
)

# Now we will merge the enrichment scores
all(rownames(gse184316_hp_ctrl_es) == rownames(gse150910_hp_ctrl_es)) # Must be TRUE

merged_score_chp_ctrl <- merge(gse184316_hp_ctrl_es,
  gse150910_hp_ctrl_es,
  by = 0,
  sort = FALSE
) %>%
  column_to_rownames(var = "Row.names")
write.csv(merged_score_chp_ctrl, "./results/xCell_results/Merged_scores_HP_Ctrl.csv",
  row.names = TRUE
)

### Reorder
merged_metadata_chp_ctrl <- merged_metadata_chp_ctrl[match(
  colnames(merged_score_chp_ctrl),
  merged_metadata_chp_ctrl$sample
), ]
write.csv(merged_metadata_chp_ctrl, "./results/xCell_results/Merged_metadata_HP_Ctrl.csv",
  row.names = TRUE
)

all(colnames(merged_score_chp_ctrl) == merged_metadata_chp_ctrl$sample) # Must be TRUE

col_data <- data.frame(
  row.names = colnames(merged_score_chp_ctrl),
  condition = factor(merged_metadata_chp_ctrl$condition,
    levels = c("Control", "HP")
  ),
  platform = factor(c(
    rep("GSE184316", ncol(gse184316_hp_ctrl_es)),
    rep("GSE150910", ncol(gse150910_hp_ctrl_es))
  ))
)

# ==============================================================================
# PCA BEFORE CORRECTION
# ==============================================================================

# ==============================================================================
# 1. SETUP & VARIANCE CALCULATION
# ==============================================================================

# Ensure PCA is run
# There cell types with enrichment scores of 0 hence it is important to remove them
# The below snippet dynamically removes those columns
# Identify rows (cell types) with zero variance
# We use apply() over rows (1) to calculate variance
# na.rm = TRUE is a safeguard in case of missing values
row_variances <- apply(merged_score_chp_ctrl, 1, var, na.rm = TRUE)

# 2. Subset the matrix to keep only rows with variance > 0
merged_score_filtered <- merged_score_chp_ctrl[row_variances > 0, ]

# 3. Verify removal (Optional)
print(paste("Removed", nrow(merged_score_chp_ctrl) - nrow(merged_score_filtered), "zero-variance features."))

# 4. Run PCA on the filtered data
pca_res <- prcomp(t(merged_score_filtered), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$platform <- col_data$platform
pca_df$condition <- col_data$condition

# Calculate Percentage Variance Explained
# summary(pca_res)$importance returns a matrix where row 2 is "Proportion of Variance"
pca_summary <- summary(pca_res)$importance
pc1_lab <- paste0("PC1 (", round(pca_summary[2, "PC1"] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(pca_summary[2, "PC2"] * 100, 1), "%)")


# ==============================================================================
# 2. PLOT 1: COLORED BY PLATFORM
# ==============================================================================

p_platform <- ggplot(pca_df, aes(x = PC1, y = PC2, color = platform)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = platform), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA: Platform Effect",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw()

# ==============================================================================
# 3. PLOT 2: COLORED BY CONDITION
# ==============================================================================

p_condition <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = condition), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA: Biological Condition",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw()

# ==============================================================================
# 4. SAVE PLOTS TO FILE
# ==============================================================================

# Option A: Save as separate files (Best for manuscript figures)
ggsave("./results/xCell_results/HP_Ctrl_score_PCA_Platform_Before.pdf",
  plot = p_platform,
  width = 16, height = 14, units = "cm", dpi = 300
)

ggsave("./results/xCell_results/HP_Ctrl_score_PCA_Condition_Before.pdf",
  plot = p_condition,
  width = 16, height = 14, units = "cm", dpi = 300
)

# Option B: Save as a combined side-by-side plot (Best for quick comparison)
# Note: we need to pass the grid object to ggsave
combined_plot <- arrangeGrob(p_platform, p_condition, ncol = 2)

ggsave("./results/xCell_results/HP_Ctrl_score_PCA_Combined_Before.pdf",
  plot = combined_plot,
  width = 32, height = 14, units = "cm", dpi = 300
)


# ==============================================================================
# COMBAT BATCH CORRECTION (Condition 1: cHP vs Control)
# ==============================================================================
# ComBat is applied directly to the merged enrichment scores (cell types x samples).
# Unlike raw counts, xCell scores are already continuous values, so no VST is needed.
# We remove platform (batch) effects while protecting the biological condition.

message("--- Running ComBat for Condition 1: cHP vs Control ---")

# 1. Prepare Inputs for ComBat
# Use the zero-variance-filtered matrix (same features used for pre-correction PCA)
score_mat <- as.matrix(merged_score_filtered)

# Define the Batch (Platform)
batch_var <- col_data$platform

# Define the Model Matrix (Biology to PROTECT)
# This tells ComBat: "Remove platform effects, but preserve differences between Control and HP"
mod_combat <- model.matrix(~condition, data = col_data)

# 2. Run ComBat
combat_scores_chp_ctrl <- ComBat(
  dat = score_mat,
  batch = batch_var,
  mod = mod_combat,
  par.prior = TRUE,
  prior.plots = FALSE
)

# 3. Export corrected enrichment scores
corrected_df <- as.data.frame(combat_scores_chp_ctrl)
write.csv(corrected_df,
  "./results/xCell_results/Corrected_scores/Corrected_scores_HP_Ctrl.csv",
  row.names = TRUE
)
message("Corrected enrichment scores saved for cHP vs Control.")

# ==============================================================================
# POST-CORRECTION PCA (Condition 1: cHP vs Control)
# ==============================================================================

pca_res_corr <- prcomp(t(combat_scores_chp_ctrl), scale. = TRUE)
pca_df_corr <- as.data.frame(pca_res_corr$x)
pca_df_corr$platform <- col_data$platform
pca_df_corr$condition <- col_data$condition

pca_summary_corr <- summary(pca_res_corr)$importance
pc1_lab_corr <- paste0("PC1 (", round(pca_summary_corr[2, "PC1"] * 100, 1), "%)")
pc2_lab_corr <- paste0("PC2 (", round(pca_summary_corr[2, "PC2"] * 100, 1), "%)")

# Plot 1: Colored by Platform (After Correction)
p_platform_corr <- ggplot(pca_df_corr, aes(x = PC1, y = PC2, color = platform)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = platform), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA After ComBat: Platform (cHP vs Ctrl)",
    x = pc1_lab_corr,
    y = pc2_lab_corr
  ) +
  theme_bw()

# Plot 2: Colored by Condition (After Correction)
p_condition_corr <- ggplot(pca_df_corr, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = condition), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA After ComBat: Condition (cHP vs Ctrl)",
    x = pc1_lab_corr,
    y = pc2_lab_corr
  ) +
  theme_bw()

# Save Post-Correction Plots
ggsave("./results/xCell_results/HP_Ctrl_score_PCA_Platform_AfterCombat.pdf",
  plot = p_platform_corr,
  width = 16, height = 14, units = "cm", dpi = 300
)

ggsave("./results/xCell_results/HP_Ctrl_score_PCA_Condition_AfterCombat.pdf",
  plot = p_condition_corr,
  width = 16, height = 14, units = "cm", dpi = 300
)

combined_plot_corr <- arrangeGrob(p_platform_corr, p_condition_corr, ncol = 2)

ggsave("./results/xCell_results/HP_Ctrl_score_PCA_Combined_AfterCombat.pdf",
  plot = combined_plot_corr,
  width = 32, height = 14, units = "cm", dpi = 300
)

message("Post-correction PCA plots saved for cHP vs Control.")

################################################################################

## Condition 2: cHP vs IPF #############################
# First we will merge metadata
gse184316_metadata_chp_ipf <- read.csv("./data/GSE184316/Sample_information_HP_IPF_GSE184316.csv",
  header = TRUE,
  row.names = 1
)

gse150910_metadata_chp_ipf <- read.csv("./data/GSE150910/Sample_information_HP_IPF_GSE150910.csv",
  header = TRUE,
  row.names = 1
)

merged_metadata_chp_ipf <- rbind(
  gse184316_metadata_chp_ipf,
  gse150910_metadata_chp_ipf
)

# Now we will merge the enrichment scores
all(rownames(gse184316_hp_ipf_es) == rownames(gse150910_hp_ipf_es)) # Must be TRUE

merged_score_chp_ipf <- merge(gse184316_hp_ipf_es,
  gse150910_hp_ipf_es,
  by = 0,
  sort = FALSE
) %>%
  column_to_rownames(var = "Row.names")
write.csv(merged_score_chp_ipf, "./results/xCell_results/Merged_scores_HP_IPF.csv",
  row.names = TRUE
)

### Reorder
merged_metadata_chp_ipf <- merged_metadata_chp_ipf[match(
  colnames(merged_score_chp_ipf),
  merged_metadata_chp_ipf$sample
), ]
write.csv(merged_metadata_chp_ipf, "./results/xCell_results/Merged_metadata_HP_IPF.csv",
  row.names = TRUE
)

all(colnames(merged_score_chp_ipf) == merged_metadata_chp_ipf$sample) # Must be TRUE

col_data <- data.frame(
  row.names = colnames(merged_score_chp_ipf),
  condition = factor(merged_metadata_chp_ipf$condition,
    levels = c("IPF", "HP")
  ),
  platform = factor(c(
    rep("GSE184316", ncol(gse184316_hp_ipf_es)),
    rep("GSE150910", ncol(gse150910_hp_ipf_es))
  ))
)

# ==============================================================================
# PCA BEFORE CORRECTION
# ==============================================================================

# ==============================================================================
# 1. SETUP & VARIANCE CALCULATION
# ==============================================================================

# Ensure PCA is run
# There cell types with enrichment scores of 0 hence it is important to remove them
# The below snippet dynamically removes those columns
# Identify rows (cell types) with zero variance
# We use apply() over rows (1) to calculate variance
# na.rm = TRUE is a safeguard in case of missing values
row_variances <- apply(merged_score_chp_ipf, 1, var, na.rm = TRUE)

# 2. Subset the matrix to keep only rows with variance > 0
merged_score_filtered <- merged_score_chp_ipf[row_variances > 0, ]

# 3. Verify removal (Optional)
print(paste("Removed", nrow(merged_score_chp_ipf) - nrow(merged_score_filtered), "zero-variance features."))

# 4. Run PCA on the filtered data
pca_res <- prcomp(t(merged_score_filtered), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$platform <- col_data$platform
pca_df$condition <- col_data$condition

# Calculate Percentage Variance Explained
# summary(pca_res)$importance returns a matrix where row 2 is "Proportion of Variance"
pca_summary <- summary(pca_res)$importance
pc1_lab <- paste0("PC1 (", round(pca_summary[2, "PC1"] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(pca_summary[2, "PC2"] * 100, 1), "%)")


# ==============================================================================
# 2. PLOT 1: COLORED BY PLATFORM
# ==============================================================================

p_platform <- ggplot(pca_df, aes(x = PC1, y = PC2, color = platform)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = platform), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA: Platform Effect",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw()

# ==============================================================================
# 3. PLOT 2: COLORED BY CONDITION
# ==============================================================================

p_condition <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = condition), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA: Biological Condition",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw()

# ==============================================================================
# 4. SAVE PLOTS TO FILE
# ==============================================================================

# Option A: Save as separate files (Best for manuscript figures)
ggsave("./results/xCell_results/HP_IPF_score_PCA_Platform_Before.pdf",
  plot = p_platform,
  width = 16, height = 14, units = "cm", dpi = 300
)

ggsave("./results/xCell_results/HP_IPF_score_PCA_Condition_Before.pdf",
  plot = p_condition,
  width = 16, height = 14, units = "cm", dpi = 300
)

# Option B: Save as a combined side-by-side plot (Best for quick comparison)
# Note: we need to pass the grid object to ggsave
combined_plot <- arrangeGrob(p_platform, p_condition, ncol = 2)

ggsave("./results/xCell_results/HP_IPF_score_PCA_Combined_Before.pdf",
  plot = combined_plot,
  width = 32, height = 14, units = "cm", dpi = 300
)


# ==============================================================================
# COMBAT BATCH CORRECTION (Condition 2: cHP vs IPF)
# ==============================================================================

message("--- Running ComBat for Condition 2: cHP vs IPF ---")

# 1. Prepare Inputs for ComBat
score_mat <- as.matrix(merged_score_filtered)
batch_var <- col_data$platform
mod_combat <- model.matrix(~condition, data = col_data)

# 2. Run ComBat
combat_scores_chp_ipf <- ComBat(
  dat = score_mat,
  batch = batch_var,
  mod = mod_combat,
  par.prior = TRUE,
  prior.plots = FALSE
)

# 3. Export corrected enrichment scores
corrected_df <- as.data.frame(combat_scores_chp_ipf)
write.csv(corrected_df,
  "./results/xCell_results/Corrected_scores/Corrected_scores_HP_IPF.csv",
  row.names = TRUE
)
message("Corrected enrichment scores saved for cHP vs IPF.")

# ==============================================================================
# POST-CORRECTION PCA (Condition 2: cHP vs IPF)
# ==============================================================================

pca_res_corr <- prcomp(t(combat_scores_chp_ipf), scale. = TRUE)
pca_df_corr <- as.data.frame(pca_res_corr$x)
pca_df_corr$platform <- col_data$platform
pca_df_corr$condition <- col_data$condition

pca_summary_corr <- summary(pca_res_corr)$importance
pc1_lab_corr <- paste0("PC1 (", round(pca_summary_corr[2, "PC1"] * 100, 1), "%)")
pc2_lab_corr <- paste0("PC2 (", round(pca_summary_corr[2, "PC2"] * 100, 1), "%)")

# Plot 1: Colored by Platform (After Correction)
p_platform_corr <- ggplot(pca_df_corr, aes(x = PC1, y = PC2, color = platform)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = platform), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA After ComBat: Platform (cHP vs IPF)",
    x = pc1_lab_corr,
    y = pc2_lab_corr
  ) +
  theme_bw()

# Plot 2: Colored by Condition (After Correction)
p_condition_corr <- ggplot(pca_df_corr, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = condition), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA After ComBat: Condition (cHP vs IPF)",
    x = pc1_lab_corr,
    y = pc2_lab_corr
  ) +
  theme_bw()

# Save Post-Correction Plots
ggsave("./results/xCell_results/HP_IPF_score_PCA_Platform_AfterCombat.pdf",
  plot = p_platform_corr,
  width = 16, height = 14, units = "cm", dpi = 300
)

ggsave("./results/xCell_results/HP_IPF_score_PCA_Condition_AfterCombat.pdf",
  plot = p_condition_corr,
  width = 16, height = 14, units = "cm", dpi = 300
)

combined_plot_corr <- arrangeGrob(p_platform_corr, p_condition_corr, ncol = 2)

ggsave("./results/xCell_results/HP_IPF_score_PCA_Combined_AfterCombat.pdf",
  plot = combined_plot_corr,
  width = 32, height = 14, units = "cm", dpi = 300
)

message("Post-correction PCA plots saved for cHP vs IPF.")


###############################################################################


## Condition 3: IPF vs Control #############################
# First we will merge metadata
gse184316_metadata_ipf_ctrl <- read.csv("./data/GSE184316/Sample_information_IPF_Control_GSE184316.csv",
  header = TRUE,
  row.names = 1
)

gse150910_metadata_ipf_ctrl <- read.csv("./data/GSE150910/Sample_information_IPF_Control_GSE150910.csv",
  header = TRUE,
  row.names = 1
)

merged_metadata_ipf_ctrl <- rbind(
  gse184316_metadata_ipf_ctrl,
  gse150910_metadata_ipf_ctrl
)

# Now we will merge the enrichment scores
all(rownames(gse184316_ipf_ctrl_es) == rownames(gse150910_ipf_ctrl_es)) # Must be TRUE

merged_score_ipf_ctrl <- merge(gse184316_ipf_ctrl_es,
  gse150910_ipf_ctrl_es,
  by = 0,
  sort = FALSE
) %>%
  column_to_rownames(var = "Row.names")
write.csv(merged_score_ipf_ctrl, "./results/xCell_results/Merged_scores_IPF_Ctrl.csv",
  row.names = TRUE
)

### Reorder
merged_metadata_ipf_ctrl <- merged_metadata_ipf_ctrl[match(
  colnames(merged_score_ipf_ctrl),
  merged_metadata_ipf_ctrl$sample
), ]
write.csv(merged_metadata_ipf_ctrl, "./results/xCell_results/Merged_metadata_IPF_Ctrl.csv",
  row.names = TRUE
)

all(colnames(merged_score_ipf_ctrl) == merged_metadata_ipf_ctrl$sample) # Must be TRUE

col_data <- data.frame(
  row.names = colnames(merged_score_ipf_ctrl),
  condition = factor(merged_metadata_ipf_ctrl$condition,
    levels = c("Control", "IPF")
  ),
  platform = factor(c(
    rep("GSE184316", ncol(gse184316_ipf_ctrl_es)),
    rep("GSE150910", ncol(gse150910_ipf_ctrl_es))
  ))
)

# ==============================================================================
# PCA BEFORE CORRECTION
# ==============================================================================

# ==============================================================================
# 1. SETUP & VARIANCE CALCULATION
# ==============================================================================

# Ensure PCA is run
# There cell types with enrichment scores of 0 hence it is important to remove them
# The below snippet dynamically removes those columns
# Identify rows (cell types) with zero variance
# We use apply() over rows (1) to calculate variance
# na.rm = TRUE is a safeguard in case of missing values
row_variances <- apply(merged_score_ipf_ctrl, 1, var, na.rm = TRUE)

# 2. Subset the matrix to keep only rows with variance > 0
merged_score_filtered <- merged_score_ipf_ctrl[row_variances > 0, ]

# 3. Verify removal (Optional)
print(paste("Removed", nrow(merged_score_ipf_ctrl) - nrow(merged_score_filtered), "zero-variance features."))

# 4. Run PCA on the filtered data
pca_res <- prcomp(t(merged_score_filtered), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$platform <- col_data$platform
pca_df$condition <- col_data$condition

# Calculate Percentage Variance Explained
# summary(pca_res)$importance returns a matrix where row 2 is "Proportion of Variance"
pca_summary <- summary(pca_res)$importance
pc1_lab <- paste0("PC1 (", round(pca_summary[2, "PC1"] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(pca_summary[2, "PC2"] * 100, 1), "%)")


# ==============================================================================
# 2. PLOT 1: COLORED BY PLATFORM
# ==============================================================================

p_platform <- ggplot(pca_df, aes(x = PC1, y = PC2, color = platform)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = platform), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA: Platform Effect",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw()

# ==============================================================================
# 3. PLOT 2: COLORED BY CONDITION
# ==============================================================================

p_condition <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = condition), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA: Biological Condition",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw()

# ==============================================================================
# 4. SAVE PLOTS TO FILE
# ==============================================================================

# Option A: Save as separate files (Best for manuscript figures)
ggsave("./results/xCell_results/IPF_Ctrl_score_PCA_Platform_Before.pdf",
  plot = p_platform,
  width = 16, height = 14, units = "cm", dpi = 300
)

ggsave("./results/xCell_results/IPF_Ctrl_score_PCA_Condition_Before.pdf",
  plot = p_condition,
  width = 16, height = 14, units = "cm", dpi = 300
)

# Option B: Save as a combined side-by-side plot (Best for quick comparison)
# Note: we need to pass the grid object to ggsave
combined_plot <- arrangeGrob(p_platform, p_condition, ncol = 2)

ggsave("./results/xCell_results/IPF_Ctrl_score_PCA_Combined_Before.pdf",
  plot = combined_plot,
  width = 32, height = 14, units = "cm", dpi = 300
)


# ==============================================================================
# COMBAT BATCH CORRECTION (Condition 3: IPF vs Control)
# ==============================================================================

message("--- Running ComBat for Condition 3: IPF vs Control ---")

# 1. Prepare Inputs for ComBat
score_mat <- as.matrix(merged_score_filtered)
batch_var <- col_data$platform
mod_combat <- model.matrix(~condition, data = col_data)

# 2. Run ComBat
combat_scores_ipf_ctrl <- ComBat(
  dat = score_mat,
  batch = batch_var,
  mod = mod_combat,
  par.prior = TRUE,
  prior.plots = FALSE
)

# 3. Export corrected enrichment scores
corrected_df <- as.data.frame(combat_scores_ipf_ctrl)
write.csv(corrected_df,
  "./results/xCell_results/Corrected_scores/Corrected_scores_IPF_Ctrl.csv",
  row.names = TRUE
)
message("Corrected enrichment scores saved for IPF vs Control.")

# ==============================================================================
# POST-CORRECTION PCA (Condition 3: IPF vs Control)
# ==============================================================================

pca_res_corr <- prcomp(t(combat_scores_ipf_ctrl), scale. = TRUE)
pca_df_corr <- as.data.frame(pca_res_corr$x)
pca_df_corr$platform <- col_data$platform
pca_df_corr$condition <- col_data$condition

pca_summary_corr <- summary(pca_res_corr)$importance
pc1_lab_corr <- paste0("PC1 (", round(pca_summary_corr[2, "PC1"] * 100, 1), "%)")
pc2_lab_corr <- paste0("PC2 (", round(pca_summary_corr[2, "PC2"] * 100, 1), "%)")

# Plot 1: Colored by Platform (After Correction)
p_platform_corr <- ggplot(pca_df_corr, aes(x = PC1, y = PC2, color = platform)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = platform), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA After ComBat: Platform (IPF vs Ctrl)",
    x = pc1_lab_corr,
    y = pc2_lab_corr
  ) +
  theme_bw()

# Plot 2: Colored by Condition (After Correction)
p_condition_corr <- ggplot(pca_df_corr, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = condition), geom = "polygon", alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "PCA After ComBat: Condition (IPF vs Ctrl)",
    x = pc1_lab_corr,
    y = pc2_lab_corr
  ) +
  theme_bw()

# Save Post-Correction Plots
ggsave("./results/xCell_results/IPF_Ctrl_score_PCA_Platform_AfterCombat.pdf",
  plot = p_platform_corr,
  width = 16, height = 14, units = "cm", dpi = 300
)

ggsave("./results/xCell_results/IPF_Ctrl_score_PCA_Condition_AfterCombat.pdf",
  plot = p_condition_corr,
  width = 16, height = 14, units = "cm", dpi = 300
)

combined_plot_corr <- arrangeGrob(p_platform_corr, p_condition_corr, ncol = 2)

ggsave("./results/xCell_results/IPF_Ctrl_score_PCA_Combined_AfterCombat.pdf",
  plot = combined_plot_corr,
  width = 32, height = 14, units = "cm", dpi = 300
)

message("Post-correction PCA plots saved for IPF vs Control.")
message("Batch correction and post-correction PCA complete for all conditions.")
