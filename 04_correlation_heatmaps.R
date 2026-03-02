# ==============================================================================
# 04_correlation_heatmaps.R
# Feature (gene-gene) and sample-sample correlation heatmaps for
# batch-corrected expression matrices from results/Corrected_matrix/
# ==============================================================================

# Load libraries
if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!require("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# ==============================================================================
# OUTPUT DIRECTORY
# ==============================================================================
out_dir <- "./results/Correlation_Heatmaps/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# COLOR PALETTE
# ==============================================================================

# Diverging blue-white-red palette for correlation values (-1 to 1)
cor_colors <- colorRampPalette(c(
    "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "white",
    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
))(100)

# Annotation colors
ann_colors <- list(
    Condition = c("Control" = "#79AC78", "HP" = "#D25353", "IPF" = "#5D688A")
)

# ==============================================================================
# CONFIGURATION: Corrected matrices + their metadata
# ==============================================================================

datasets <- list(
    list(
        name      = "cHP_vs_Ctrl",
        expr_file = "./results/Corrected_matrix/Merged_cHPvsCtrl_BatchCorrected_Expression.csv",
        meta_file = "./results/xCell_results/Merged_metadata_HP_Ctrl.csv"
    ),
    list(
        name      = "cHP_vs_IPF",
        expr_file = "./results/Corrected_matrix/Merged_cHPvsIPF_BatchCorrected_Expression.csv",
        meta_file = "./results/xCell_results/Merged_metadata_HP_IPF.csv"
    ),
    list(
        name      = "IPF_vs_Ctrl",
        expr_file = "./results/Corrected_matrix/Merged_IPFvsCtrl_BatchCorrected_Expression.csv",
        meta_file = "./results/xCell_results/Merged_metadata_IPF_Ctrl.csv"
    )
)

# ==============================================================================
# MAIN LOOP
# ==============================================================================

for (ds in datasets) {
    message(paste0("\n===== Processing: ", ds$name, " ====="))

    # --------------------------------------------------------------------------
    # 1. LOAD EXPRESSION MATRIX (genes x samples)
    # --------------------------------------------------------------------------
    if (!file.exists(ds$expr_file)) {
        message(paste("  File not found:", ds$expr_file, "- skipping."))
        next
    }

    expr <- read.csv(ds$expr_file, row.names = 1, check.names = FALSE)
    expr <- as.matrix(expr)
    message(paste("  Loaded:", nrow(expr), "genes x", ncol(expr), "samples"))

    # --------------------------------------------------------------------------
    # 2. BUILD SAMPLE ANNOTATION (condition only)
    # --------------------------------------------------------------------------
    meta <- read.csv(ds$meta_file)

    # Align metadata to expression columns
    meta <- meta[match(colnames(expr), meta$sample), ]

    # Create annotation data frame for pheatmap (Condition only)
    sample_ann <- data.frame(
        Condition = meta$condition,
        row.names = colnames(expr)
    )

    # --------------------------------------------------------------------------
    # 3. SAMPLE-SAMPLE CORRELATION HEATMAP
    # --------------------------------------------------------------------------
    message("  Computing sample-sample correlation...")

    sample_cor <- cor(expr, method = "pearson")

    # Plot
    p_sample <- pheatmap(
        sample_cor,
        color              = cor_colors,
        breaks             = seq(-1, 1, length.out = 101),
        clustering_method  = "ward.D2",
        annotation_col     = sample_ann,
        annotation_row     = sample_ann,
        annotation_colors  = ann_colors,
        show_rownames      = FALSE,
        show_colnames      = FALSE,
        fontsize           = 8,
        border_color       = NA,
        main               = paste0("Sample Correlation \u2014 ", ds$name, " (Batch Corrected)"),
        silent             = TRUE
    )

    # Determine dynamic size (capped to stay within pixel limits at 300 dpi)
    n_samples <- ncol(expr)
    sample_size <- min(40, max(8, n_samples * 0.08 + 3))

    # Save PNG (300 dpi) using direct device call
    png_file <- paste0(out_dir, "Sample_Correlation_", ds$name, ".png")
    png(png_file, width = sample_size, height = sample_size, units = "in", res = 300)
    grid::grid.newpage()
    grid::grid.draw(p_sample$gtable)
    dev.off()
    message(paste("  Saved:", png_file))

    # Save PDF (vector — no pixel concern)
    pdf_file <- paste0(out_dir, "Sample_Correlation_", ds$name, ".pdf")
    pdf(pdf_file, width = sample_size, height = sample_size)
    grid::grid.newpage()
    grid::grid.draw(p_sample$gtable)
    dev.off()
    message(paste("  Saved:", pdf_file))

    # --------------------------------------------------------------------------
    # 4. FEATURE (GENE-GENE) CORRELATION HEATMAP
    # --------------------------------------------------------------------------
    # Top 100 most variable genes are selected using row-wise variance:
    #   var(gene_i) = variance of that gene's expression across ALL samples.
    # Genes are ranked by variance (descending) and the top 100 are kept.
    # High-variance genes drive biological separation; low-variance genes
    # add noise and make the heatmap unreadable.
    message("  Computing feature (gene-gene) correlation...")

    gene_vars <- apply(expr, 1, var)
    max_features <- 100

    if (nrow(expr) > max_features) {
        top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:max_features]
        expr_sub <- expr[top_genes, ]
        feat_title <- paste0(
            "Top ", max_features, " Variable Genes \u2014 Feature Correlation \u2014 ",
            ds$name, " (Batch Corrected)"
        )
        message(paste("  Subset to top", max_features, "variable genes."))
    } else {
        expr_sub <- expr
        feat_title <- paste0("Feature Correlation \u2014 ", ds$name, " (Batch Corrected)")
    }

    feature_cor <- cor(t(expr_sub), method = "pearson")

    # Plot
    p_feature <- pheatmap(
        feature_cor,
        color              = cor_colors,
        breaks             = seq(-1, 1, length.out = 101),
        clustering_method  = "ward.D2",
        show_rownames      = TRUE,
        show_colnames      = FALSE,
        fontsize           = 8,
        fontsize_row       = 4,
        border_color       = NA,
        main               = feat_title,
        silent             = TRUE
    )

    # Determine dynamic size
    n_features <- nrow(expr_sub)
    feature_size <- min(25, max(10, n_features * 0.12 + 3))

    # Save PNG (300 dpi) using direct device call
    png_file <- paste0(out_dir, "Feature_Correlation_", ds$name, ".png")
    png(png_file, width = feature_size, height = feature_size, units = "in", res = 300)
    grid::grid.newpage()
    grid::grid.draw(p_feature$gtable)
    dev.off()
    message(paste("  Saved:", png_file))

    # Save PDF (vector)
    pdf_file <- paste0(out_dir, "Feature_Correlation_", ds$name, ".pdf")
    pdf(pdf_file, width = feature_size, height = feature_size)
    grid::grid.newpage()
    grid::grid.draw(p_feature$gtable)
    dev.off()
    message(paste("  Saved:", pdf_file))
}

message("\n===== All correlation heatmaps complete. =====")
