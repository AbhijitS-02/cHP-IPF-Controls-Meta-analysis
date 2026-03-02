# Load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("sva")

library(tidyverse)
library(DESeq2)
library(sva)
library(ggplot2)

out_dir <- "./results/PCA_plots/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("./results/Corrected_matrix/", showWarnings = FALSE, recursive = TRUE)
dir.create("./results/Raw_merged_data", showWarnings = FALSE, recursive = TRUE)
# Read the FINM metaDEGs that are selected according to the following criteria:
# 1) mean absolute logFC > 1
# 2) meta p-adjusted < 0.05
# 3) must be expressed in all the 3 datasets
# 4) must have the same direction (+ or -) of expression across all the datasets

chp_ctrl <- read.csv("cHP_Ctrl_Common_genes_same_direction_expr.csv", 
                     header = TRUE)

chp_ipf <- read.csv("cHP_IPF_Common_genes_same_direction_expr.csv",
                    header = TRUE)

ipf_ctrl <- read.csv("IPF_Ctrl_Common_genes_same_direction_expr.csv",
                     header = TRUE)

# Find union of genes (unique and intersecting among the datasets - has no duplicates)
union_genes <- Reduce(union, list(chp_ctrl$Gene_symbol,
                                  chp_ipf$Gene_symbol,
                                  ipf_ctrl$Gene_symbol))

# Merge datasets for 3 conditions:
#################################
# Condition 1) cHP vs Controls
GSE184316_chp_ctrl <- read.csv("data/GSE184316/GSE184316_raw_counts_HP_Control.csv",
                               header = TRUE,
                               row.names = 1)

GSE150910_chp_ctrl <- read.csv("data/GSE150910/GSE150910_raw_counts_HP_Control.csv",
                               header = TRUE,
                               row.names = 1)

## Subset
common <- intersect(rownames(GSE184316_chp_ctrl), union_genes)
GSE184316_chp_ctrl <- GSE184316_chp_ctrl[common, ]

common <- intersect(rownames(GSE150910_chp_ctrl), union_genes)
GSE150910_chp_ctrl <- GSE150910_chp_ctrl[common,]


# Condition 2) cHP vs IPF
GSE184316_chp_ipf <- read.csv("data/GSE184316/GSE184316_raw_counts_HP_IPF.csv",
                              header = TRUE,
                              row.names = 1)

GSE150910_chp_ipf <- read.csv("data/GSE150910/GSE150910_raw_counts_HP_IPF.csv",
                              header = TRUE,
                              row.names = 1)

## Subset
common <- intersect(rownames(GSE184316_chp_ipf), union_genes)
GSE184316_chp_ipf <- GSE184316_chp_ipf[common, ]

common <- intersect(rownames(GSE150910_chp_ipf), union_genes)
GSE150910_chp_ipf <- GSE150910_chp_ipf[common,]



# Condtion 3) IPF vs Control
GSE184316_ipf_ctrl <- read.csv("data/GSE184316/GSE184316_raw_counts_IPF_Control.csv",
                               header = TRUE,
                               row.names = 1)

GSE150910_ipf_ctrl <- read.csv("data/GSE150910/GSE150910_raw_counts_IPF_Control.csv",
                               header = TRUE,
                               row.names = 1)

## Subset
common <- intersect(rownames(GSE184316_ipf_ctrl), union_genes)
GSE184316_ipf_ctrl <- GSE184316_ipf_ctrl[common, ]

common <- intersect(rownames(GSE150910_ipf_ctrl), union_genes)
GSE150910_ipf_ctrl <- GSE150910_ipf_ctrl[common,]

# Now that subsetting is done
# We will find if there are any NA values and remove them
# Total count of NAs across the entire dataframe
print(paste("Total NA values in GSE184316 cHP vs Control:", sum(is.na(GSE184316_chp_ctrl))))
print(paste("Total NA values in GSE184316 cHP vs IPF:", sum(is.na(GSE184316_chp_ipf))))
print(paste("Total NA values in GSE184316 IPF vs Control:", sum(is.na(GSE184316_ipf_ctrl))))

print(paste("Total NA values in GSE150910 cHP vs Control:", sum(is.na(GSE150910_chp_ctrl))))
print(paste("Total NA values in GSE150910 cHP vs IPF:", sum(is.na(GSE150910_chp_ipf))))
print(paste("Total NA values in GSE150910 IPF vs Control:", sum(is.na(GSE150910_ipf_ctrl))))


# Now we will pairwise merge each condition to create merged dataset
# eg: Pairwise = GSE184316_chp_ipf + GSE150910_chp_ipf

# Condition 1: cHP vs Control
merged_chp_ctrl <- merge(GSE184316_chp_ctrl,
                        GSE150910_chp_ctrl,
                        by = 0,
                        sort = FALSE) %>%
  column_to_rownames(var = "Row.names")
write.csv(merged_chp_ctrl, 
          "./results/Raw_merged_data/Merged_raw_HP_Ctrl.csv",
          row.names = TRUE)

# Condition 2: cHP vs IPF
merged_chp_ipf <- merge(GSE184316_chp_ipf,
                        GSE150910_chp_ipf,
                        by = 0,
                        sort = FALSE) %>%
  column_to_rownames(var = "Row.names")
write.csv(merged_chp_ipf,
          "./results/Raw_merged_data/Merged_raw_HP_IPF.csv",
          row.names = TRUE)


# Condition 3: IPF vs Control
merged_ipf_ctrl <- merge(GSE184316_ipf_ctrl,
                        GSE150910_ipf_ctrl,
                        by = 0,
                        sort = FALSE) %>%
  column_to_rownames(var = "Row.names")
write.csv(merged_ipf_ctrl,
          "./results/Raw_merged_data/Merged_raw_IPF_Ctrl.csv",
          row.names = TRUE)

# Check if NA present
print(paste("Total NA values in Merged cHP vs Control:", sum(is.na(merged_chp_ctrl))))
print(paste("Total NA values in Merged cHP vs IPF:", sum(is.na(merged_chp_ipf))))
print(paste("Total NA values in Merged IPF vs Control:", sum(is.na(merged_ipf_ctrl))))

# Since no NA values present, we move on.....
###########################################################################

# Plot PCA for the 3 conditions
my_colors <- c("Control"   =   "#79AC78", 
               "HP"        =   "#D25353", 
               "IPF"       =   "#5D688A",
               "GSE184316" =   "#EF9C66",
               "GSE150910" =   "#9ECAD6")


## Condition 1: cHP vs Control
message("--- Running PCA for Condition 1: cHP vs Control ---")
gse184316_metadata_chp_ctrl <- read.csv("./data/GSE184316/Sample_information_HP_Control_GSE184316.csv",
                                       header = TRUE,
                                       row.names = 1)

gse150910_metadata_chp_ctrl <- read.csv("./data/GSE150910/Sample_information_HP_Control_GSE150910.csv",
                                       header = TRUE,
                                       row.names = 1)

merged_metadata_chp_ctrl <- rbind(gse184316_metadata_chp_ctrl,
                                  gse150910_metadata_chp_ctrl)

### Reorder
merged_metadata_chp_ctrl <- merged_metadata_chp_ctrl[match(colnames(merged_chp_ctrl), 
                                                           merged_metadata_chp_ctrl$sample),]
all(colnames(merged_chp_ctrl) == merged_metadata_chp_ctrl$sample) # Must be TRUE


# Create colData
# 'condition' is the biological variable, 'batch' is the dataset origin
col_data <- data.frame(
  row.names = colnames(merged_chp_ctrl),
  condition = factor(merged_metadata_chp_ctrl$condition, 
                     levels = c("Control", "HP")), 
  platform = factor(c(rep("GSE184316", ncol(GSE184316_chp_ctrl)), 
                   rep("GSE150910", ncol(GSE150910_chp_ctrl))))
)

# 3. Create DDS Object
dds <- DESeqDataSetFromMatrix(countData = merged_chp_ctrl,
                              colData = col_data,
                              design = ~ platform + condition)

# 4. VST Transform (blind=TRUE is optimal for QC/PCA)
vsd <- vst(dds, blind = TRUE)

# 5. Plot PCA
# Color by batch to visualize the specific effect you are looking for
p1 <- plotPCA(vsd, intgroup = c("platform")) + 
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = platform), alpha = 0.25) +
  ggtitle("PCA by Platform") +
  theme_bw() +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

p2 <- plotPCA(vsd, intgroup = c("condition")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = condition), alpha = 0.25) +
  ggtitle("PCA by Condition") +
  theme_bw() +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

# --- A) Save Separate PNGs (300 dpi) ---
ggsave(filename = paste0(out_dir, "Merged_HPvsCtrl_PCA_Platform.png"), plot = p1, 
       width = 8, height = 8, dpi = 300)
ggsave(filename = paste0(out_dir, "Merged_HPvsCtrl_PCA_Condition.png"), plot = p2, 
       width = 8, height = 8, dpi = 300)

# --- B) Save Separate PDFs ---
ggsave(filename = paste0(out_dir, "Merged_HPvsCtrl_PCA_Platform.pdf"), plot = p1, 
       width = 8, height = 8)
ggsave(filename = paste0(out_dir, "Merged_HPvsCtrl_PCA_Condition.pdf"), plot = p2, 
       width = 8, height = 8)

message("All plots saved successfully to ", out_dir)


# Since we observed platform effect, i.e., samples are clustered according to
# their dataset/platform (Illumina vs Ion Torrent), we must remove the technical
# differences and preserve the biological variation
message("--- Running ComBat for Condition 1: cHP vs Control ---")

# 1. Prepare Inputs for ComBat
# Extract VST-transformed matrix (Genes x Samples)
vst_mat <- assay(vsd) 

# Define the Batch (Platform)
batch_var <- colData(dds)$platform

# Define the Model Matrix (Biology to PROTECT)
# This tells ComBat: "Remove platform effects, but preserve differences between Control and HP"
mod_combat <- model.matrix(~ condition, data = colData(dds))

# 2. Run ComBat
combat_data <- ComBat(dat = vst_mat, 
                      batch = batch_var, 
                      mod = mod_combat, 
                      par.prior = TRUE, 
                      prior.plots = FALSE)

# 3. Overwrite the VST object with Corrected Data
# This allows us to use plotPCA() on the corrected data easily
vsd_corrected <- vsd
assay(vsd_corrected) <- combat_data

# 4. Plot PCA (After Correction)
p1_corr <- plotPCA(vsd_corrected, intgroup = c("platform")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = platform), alpha = 0.25) +
  ggtitle("PCA After ComBat: Platform (cHP vs Ctrl)") +
  theme_bw() + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

p2_corr <- plotPCA(vsd_corrected, intgroup = c("condition")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = condition), alpha = 0.25) +
  ggtitle("PCA After ComBat: Condition (cHP vs Ctrl)") +
  theme_bw() + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

# 5. Save Plots
ggsave(filename = paste0(out_dir, "Merged_HPvsCtrl_PCA_Platform_AfterCombat.png"), 
       plot = p1_corr, width = 8, height = 8, dpi = 300)
ggsave(filename = paste0(out_dir, "Merged_HPvsCtrl_PCA_Condition_AfterCombat.png"),
       plot = p2_corr, width = 8, height = 8, dpi = 300)

ggsave(filename = paste0(out_dir, "Merged_HPvsCtrl_PCA_Platform_AfterCombat.pdf"),
       plot = p1_corr, width = 8, height = 8)
ggsave(filename = paste0(out_dir, "Merged_HPvsCtrl_PCA_Condition_AfterCombat.pdf"),
       plot = p2_corr, width = 8, height = 8)

# Export corrected matrix
corrected_df <- as.data.frame(combat_data)
output_csv <- paste0("./results/Corrected_matrix/Merged_cHPvsCtrl_BatchCorrected_Expression.csv")
write.csv(corrected_df, file = output_csv, row.names = TRUE)
message("Batch corrected expression matrix saved to: ", output_csv)

###################################################################

# Condition 2: cHP vs IPF
message("--- Running PCA for Condition 2: cHP vs IPF ---")
gse184316_metadata_chp_ipf <- read.csv("./data/GSE184316/Sample_information_HP_IPF_GSE184316.csv",
                                        header = TRUE,
                                        row.names = 1)

gse150910_metadata_chp_ipf <- read.csv("./data/GSE150910/Sample_information_HP_IPF_GSE150910.csv",
                                        header = TRUE,
                                        row.names = 1)

merged_metadata_chp_ipf <- rbind(gse184316_metadata_chp_ipf,
                                  gse150910_metadata_chp_ipf)

### Reorder
merged_metadata_chp_ipf <- merged_metadata_chp_ipf[match(colnames(merged_chp_ipf), 
                                                           merged_metadata_chp_ipf$sample),]
all(colnames(merged_chp_ipf) == merged_metadata_chp_ipf$sample) # Must be TRUE


# Create colData
# 'condition' is the biological variable, 'batch' is the dataset origin
col_data <- data.frame(
  row.names = colnames(merged_chp_ipf),
  condition = factor(merged_metadata_chp_ipf$condition, 
                     levels = c("IPF", "HP")), 
  platform = factor(c(rep("GSE184316", ncol(GSE184316_chp_ipf)), 
                      rep("GSE150910", ncol(GSE150910_chp_ipf))))
)

# 3. Create DDS Object
dds <- DESeqDataSetFromMatrix(countData = merged_chp_ipf,
                              colData = col_data,
                              design = ~ platform + condition)

# 4. VST Transform (blind=TRUE is optimal for QC/PCA)
vsd <- vst(dds, blind = TRUE)

# 5. Plot PCA
# Color by batch to visualize the specific effect you are looking for
p1 <- plotPCA(vsd, intgroup = c("platform")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = platform), alpha = 0.25) +
  ggtitle("PCA by Platform") +
  theme_bw() +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

p2 <- plotPCA(vsd, intgroup = c("condition")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = condition), alpha = 0.25) +
  ggtitle("PCA by Condition") +
  theme_bw() +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill


# --- A) Save Separate PNGs (300 dpi) ---
ggsave(filename = paste0(out_dir, "Merged_HPvsIPF_PCA_Platform.png"), plot = p1, 
       width = 8, height = 8, dpi = 300)
ggsave(filename = paste0(out_dir, "Merged_HPvsIPF_PCA_Condition.png"), plot = p2, 
       width = 8, height = 8, dpi = 300)

# --- B) Save Separate PDFs ---
ggsave(filename = paste0(out_dir, "Merged_HPvsIPF_PCA_Platform.pdf"), plot = p1, 
       width = 8, height = 8)
ggsave(filename = paste0(out_dir, "Merged_HPvsIPF_PCA_Condition.pdf"), plot = p2, 
       width = 8, height = 8)

message("All plots saved successfully to ", out_dir)


# Since we observed platform effect, i.e., samples are clustered according to
# their dataset/platform (Illumina vs Ion Torrent), we must remove the technical
# differences and preserve the biological variation
message("--- Running ComBat for Condition 2: cHP vs IPF ---")

# 1. Prepare Inputs
vst_mat <- assay(vsd) 
batch_var <- colData(dds)$platform
mod_combat <- model.matrix(~ condition, data = colData(dds))

# 2. Run ComBat
combat_data <- ComBat(dat = vst_mat, 
                      batch = batch_var, 
                      mod = mod_combat, 
                      par.prior = TRUE, 
                      prior.plots = FALSE)

# 3. Update Object
vsd_corrected <- vsd
assay(vsd_corrected) <- combat_data

# 4. Plot
p1_corr <- plotPCA(vsd_corrected, intgroup = c("platform")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = platform), alpha = 0.25) +
  ggtitle("PCA After ComBat: Platform (cHP vs IPF)") +
  theme_bw() + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

p2_corr <- plotPCA(vsd_corrected, intgroup = c("condition")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = condition), alpha = 0.25) +
  ggtitle("PCA After ComBat: Condition (cHP vs IPF)") +
  theme_bw() + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

# 5. Save
ggsave(filename = paste0(out_dir, "Merged_HPvsIPF_PCA_Platform_AfterCombat.png"),
       plot = p1_corr, width = 8, height = 8, dpi = 300)
ggsave(filename = paste0(out_dir, "Merged_HPvsIPF_PCA_Condition_AfterCombat.png"),
       plot = p2_corr, width = 8, height = 8, dpi = 300)

ggsave(filename = paste0(out_dir, "Merged_HPvsIPF_PCA_Platform_AfterCombat.pdf"),
       plot = p1_corr, width = 8, height = 8)
ggsave(filename = paste0(out_dir, "Merged_HPvsIPF_PCA_Condition_AfterCombat.pdf"),
       plot = p2_corr, width = 8, height = 8)

# Export corrected matrix
corrected_df <- as.data.frame(combat_data)
output_csv <- paste0("./results/Corrected_matrix/Merged_cHPvsIPF_BatchCorrected_Expression.csv")
write.csv(corrected_df, file = output_csv, row.names = TRUE)
message("Batch corrected expression matrix saved to: ", output_csv)

###################################################################


# Condition 3: IPF vs Control
message("--- Running PCA for Condition 3: IPF vs Control ---")
gse184316_metadata_ipf_ctrl <- read.csv("./data/GSE184316/Sample_information_IPF_Control_GSE184316.csv",
                                       header = TRUE,
                                       row.names = 1)

gse150910_metadata_ipf_ctrl <- read.csv("./data/GSE150910/Sample_information_IPF_Control_GSE150910.csv",
                                       header = TRUE,
                                       row.names = 1)

merged_metadata_ipf_ctrl <- rbind(gse184316_metadata_ipf_ctrl,
                                 gse150910_metadata_ipf_ctrl)

### Reorder
merged_metadata_ipf_ctrl <- merged_metadata_ipf_ctrl[match(colnames(merged_ipf_ctrl), 
                                                         merged_metadata_ipf_ctrl$sample),]
all(colnames(merged_ipf_ctrl) == merged_metadata_ipf_ctrl$sample) # Must be TRUE


# Create colData
# 'condition' is the biological variable, 'batch' is the dataset origin
col_data <- data.frame(
  row.names = colnames(merged_ipf_ctrl),
  condition = factor(merged_metadata_ipf_ctrl$condition, 
                     levels = c("Control", "IPF")), 
  platform = factor(c(rep("GSE184316", ncol(GSE184316_ipf_ctrl)), 
                      rep("GSE150910", ncol(GSE150910_ipf_ctrl))))
)

# 3. Create DDS Object
dds <- DESeqDataSetFromMatrix(countData = merged_ipf_ctrl,
                              colData = col_data,
                              design = ~ platform + condition)

# 4. VST Transform (blind=TRUE is optimal for QC/PCA)
vsd <- vst(dds, blind = TRUE)

# 5. Plot PCA
# Color by batch to visualize the specific effect you are looking for
p1 <- plotPCA(vsd, intgroup = c("platform")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = platform), alpha = 0.25) +
  ggtitle("PCA by Platform") +
  theme_bw() + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

p2 <- plotPCA(vsd, intgroup = c("condition")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = condition), alpha = 0.25) +
  ggtitle("PCA by Condition") +
  theme_bw() +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill


# --- A) Save Separate PNGs (300 dpi) ---
ggsave(filename = paste0(out_dir, "Merged_IPFvsCtrl_PCA_Platform.png"), plot = p1, 
       width = 8, height = 8, dpi = 300)
ggsave(filename = paste0(out_dir, "Merged_IPFvsCtrl_PCA_Condition.png"), plot = p2, 
       width = 8, height = 8, dpi = 300)

# --- B) Save Separate PDFs ---
ggsave(filename = paste0(out_dir, "Merged_IPFvsCtrl_PCA_Platform.pdf"), plot = p1, 
       width = 8, height = 8)
ggsave(filename = paste0(out_dir, "Merged_IPFvsCtrl_PCA_Condition.pdf"), plot = p2, 
       width = 8, height = 8)

message("All plots saved successfully to ", out_dir)


# Since we observed platform effect, i.e., samples are clustered according to
# their dataset/platform (Illumina vs Ion Torrent), we must remove the technical
# differences and preserve the biological variation
message("--- Running ComBat for Condition 3: IPF vs Control ---")

# 1. Prepare Inputs
vst_mat <- assay(vsd) 
batch_var <- colData(dds)$platform
mod_combat <- model.matrix(~ condition, data = colData(dds))

# 2. Run ComBat
combat_data <- ComBat(dat = vst_mat, 
                      batch = batch_var, 
                      mod = mod_combat, 
                      par.prior = TRUE, 
                      prior.plots = FALSE)

# 3. Update Object
vsd_corrected <- vsd
assay(vsd_corrected) <- combat_data

# 4. Plot
p1_corr <- plotPCA(vsd_corrected, intgroup = c("platform")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = platform), alpha = 0.25) +
  ggtitle("PCA After ComBat: Platform (IPF vs Control)") +
  theme_bw() + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

p2_corr <- plotPCA(vsd_corrected, intgroup = c("condition")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = condition), alpha = 0.25) +
  ggtitle("PCA After ComBat: Condition (IPF vs Control)") +
  theme_bw() + 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) # Added to color the ellipse fill

# 5. Save
ggsave(filename = paste0(out_dir, "Merged_IPFvsCtrl_PCA_Platform_AfterCombat.png"),
       plot = p1_corr, width = 8, height = 8, dpi = 300)
ggsave(filename = paste0(out_dir, "Merged_IPFvsCtrl_PCA_Condition_AfterCombat.png"),
       plot = p2_corr, width = 8, height = 8, dpi = 300)

ggsave(filename = paste0(out_dir, "Merged_IPFvsCtrl_PCA_Platform_AfterCombat.pdf"),
       plot = p1_corr, width = 8, height = 8)
ggsave(filename = paste0(out_dir, "Merged_IPFvsCtrl_PCA_Condition_AfterCombat.pdf"),
       plot = p2_corr, width = 8, height = 8)

# Export corrected matrix
corrected_df <- as.data.frame(combat_data)
output_csv <- paste0("./results/Corrected_matrix/Merged_IPFvsCtrl_BatchCorrected_Expression.csv")
write.csv(corrected_df, file = output_csv, row.names = TRUE)
message("Batch corrected expression matrix saved to: ", output_csv)

message("Batch correction and plotting complete for all conditions.")

###############################################################################
# Move to 02_tpm.R
