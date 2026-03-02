################################################################################
#  08_independent_AUC_ROC.R
#
#  PURPOSE:
#    For each of the 3 conditions (cHP_Ctrl, cHP_IPF, IPF_Ctrl), take the
#    top genes (from boxplot statistics) and their combinations (from LASSO/GLM
#    results), then compute AUC-ROC on three independent raw-count datasets:
#
#      1. GSE184316          – full condition raw counts
#      2. GSE150910 (Biopsy) – biopsy-only subset of raw counts
#      3. GSE150910 (Explant)– explant-only subset of raw counts
#
#    • Single genes  : direct ROC (pROC::roc on expression vector vs label)
#    • Combinations   : logistic regression (glm) → predicted probabilities → ROC
#
#    Results are saved as separate CSVs under ./results/independent_AUC_ROC/.
#
################################################################################

# ── 0. Libraries ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(pROC)
})

# ── 1. Setup ─────────────────────────────────────────────────────────────────
cat("============================================================\n")
cat("  08_independent_AUC_ROC.R – Independent AUC-ROC Validation\n")
cat("============================================================\n\n")

# Define base directory (script's working directory)
base_dir <- "."

# Output directory
out_dir <- file.path(base_dir, "results", "independent_AUC_ROC")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── 2. Condition Mapping ─────────────────────────────────────────────────────
# Each condition has:
#   - A short label used in filenames
#   - Corresponding GSE184316 raw counts + sample info
#   - Corresponding GSE150910 raw counts + biopsy/explant sample info
#   - The condition column name and how groups are encoded
#
# GSE184316: sample info has columns: sample, condition
#   HP_Control : condition = "HP" (positive=1) vs "Control" (reference=0)
#   HP_IPF     : condition = "HP" (positive=1) vs "IPF" (reference=0)
#   IPF_Control: condition = "IPF" (positive=1) vs "Control" (reference=0)
#
# GSE150910 biopsy/explant: columns: rowname, title, Group, Sample type, Group_Encoded
#   Group_Encoded already binary (1 = positive, 0 = reference)

conditions <- list(
    list(
        label              = "cHP_Ctrl",
        gse184316_counts   = file.path(base_dir, "data", "GSE184316", "GSE184316_raw_counts_HP_Control.csv"),
        gse184316_info     = file.path(base_dir, "data", "GSE184316", "Sample_information_HP_Control_GSE184316.csv"),
        gse184316_pos      = "HP", # positive class label in 'condition' column
        gse150910_counts   = file.path(base_dir, "data", "GSE150910", "GSE150910_raw_counts_HP_Control.csv"),
        gse150910_biopsy   = file.path(base_dir, "data", "GSE150910", "Sample_info_biopsy_HP_Ctrl.csv"),
        gse150910_explant  = file.path(base_dir, "data", "GSE150910", "Sample_info_explant_HP_Ctrl.csv"),
        gene_stats_file    = file.path(base_dir, "results", "boxplots", "cHP_Ctrl", "cHP_Ctrl_gene_statistics.csv"),
        glm_combos_file    = file.path(base_dir, "results", "lasso", "cHP_Ctrl", "GLM_combinations_AUC_results.csv")
    ),
    list(
        label              = "cHP_IPF",
        gse184316_counts   = file.path(base_dir, "data", "GSE184316", "GSE184316_raw_counts_HP_IPF.csv"),
        gse184316_info     = file.path(base_dir, "data", "GSE184316", "Sample_information_HP_IPF_GSE184316.csv"),
        gse184316_pos      = "HP",
        gse150910_counts   = file.path(base_dir, "data", "GSE150910", "GSE150910_raw_counts_HP_IPF.csv"),
        gse150910_biopsy   = file.path(base_dir, "data", "GSE150910", "Sample_info_biopsy_HP_IPF.csv"),
        gse150910_explant  = file.path(base_dir, "data", "GSE150910", "Sample_info_explant_HP_IPF.csv"),
        gene_stats_file    = file.path(base_dir, "results", "boxplots", "cHP_IPF", "cHP_IPF_gene_statistics.csv"),
        glm_combos_file    = file.path(base_dir, "results", "lasso", "cHP_IPF", "GLM_combinations_AUC_results.csv")
    ),
    list(
        label              = "IPF_Ctrl",
        gse184316_counts   = file.path(base_dir, "data", "GSE184316", "GSE184316_raw_counts_IPF_Control.csv"),
        gse184316_info     = file.path(base_dir, "data", "GSE184316", "Sample_information_IPF_Control_GSE184316.csv"),
        gse184316_pos      = "IPF",
        gse150910_counts   = file.path(base_dir, "data", "GSE150910", "GSE150910_raw_counts_IPF_Control.csv"),
        gse150910_biopsy   = file.path(base_dir, "data", "GSE150910", "Sample_info_biopsy_IPF_Ctrl.csv"),
        gse150910_explant  = file.path(base_dir, "data", "GSE150910", "Sample_info_explant_IPF_Ctrl.csv"),
        gene_stats_file    = file.path(base_dir, "results", "boxplots", "IPF_Ctrl", "IPF_Ctrl_gene_statistics.csv"),
        glm_combos_file    = file.path(base_dir, "results", "lasso", "IPF_Ctrl", "GLM_combinations_AUC_results.csv")
    )
)

################################################################################
# ── 3. Helper Functions ──────────────────────────────────────────────────────
################################################################################

#' Prepare a dataset from GSE184316 raw counts
#' @param counts_file  Path to CSV with genes as rows, samples as columns.
#'                     First column is gene symbol / row name.
#' @param info_file    Path to sample information CSV with columns: sample, condition
#' @param pos_label    The condition label to treat as the positive class (1).
#' @return A list with:  mat = numeric matrix (samples x genes),
#'                       labels = binary vector (1=positive, 0=reference)
prepare_gse184316 <- function(counts_file, info_file, pos_label) {
    cat("    Reading GSE184316 counts:", counts_file, "\n")
    counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)

    cat("    Reading GSE184316 sample info:", info_file, "\n")
    info <- read.csv(info_file, check.names = FALSE)

    # The sample column may be quoted or not; extract it
    sample_col <- info$sample
    cond_col <- info$condition

    # Map samples to binary: positive class = 1, reference = 0
    labels <- ifelse(cond_col == pos_label, 1, 0)
    names(labels) <- sample_col

    # Subset counts to only samples in the sample info (and in the same order)
    shared_samples <- intersect(sample_col, colnames(counts))
    cat("      Shared samples (GSE184316):", length(shared_samples), "\n")

    mat <- t(counts[, shared_samples, drop = FALSE]) # samples x genes
    labels <- labels[shared_samples]

    list(mat = mat, labels = labels)
}

#' Prepare a biopsy or explant subset from GSE150910 raw counts
#' @param counts_file    Path to GSE150910 raw counts CSV (gene=symbol col, then sample cols)
#' @param sample_file    Path to biopsy or explant sample info CSV
#' @return A list with:  mat = numeric matrix (samples x genes),
#'                       labels = binary vector (Group_Encoded)
prepare_gse150910_subset <- function(counts_file, sample_file) {
    cat("    Reading GSE150910 counts:", counts_file, "\n")
    counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)

    cat("    Reading sample info:", sample_file, "\n")
    info <- read.csv(sample_file, row.names = 1, check.names = FALSE)

    # Sample IDs are the rownames of the info file
    sample_ids <- rownames(info)
    labels <- info$Group_Encoded
    names(labels) <- sample_ids

    # Subset counts to these samples
    shared_samples <- intersect(sample_ids, colnames(counts))
    cat("      Shared samples (GSE150910 subset):", length(shared_samples), "\n")

    mat <- t(counts[, shared_samples, drop = FALSE]) # samples x genes
    labels <- labels[shared_samples]

    list(mat = mat, labels = labels)
}

#' Compute single-gene AUC using pROC::roc
#' @param gene_name  Gene name (must exist as a column in mat)
#' @param mat        Numeric matrix (samples x genes)
#' @param labels     Binary label vector (same order as rows of mat)
#' @return AUC value (numeric) or NA if gene not found / error
compute_single_gene_auc <- function(gene_name, mat, labels) {
    if (!(gene_name %in% colnames(mat))) {
        return(NA)
    }
    expr <- mat[, gene_name]
    tryCatch(
        {
            roc_obj <- roc(labels, expr, quiet = TRUE, direction = "auto")
            as.numeric(auc(roc_obj))
        },
        error = function(e) {
            NA
        }
    )
}

#' Compute combination-gene AUC using GLM → pROC
#' @param gene_names Character vector of gene names in the combination
#' @param mat        Numeric matrix (samples x genes)
#' @param labels     Binary label vector
#' @return AUC value (numeric) or NA if any gene missing / error
compute_combo_auc <- function(gene_names, mat, labels) {
    # Check all genes are present
    missing <- gene_names[!(gene_names %in% colnames(mat))]
    if (length(missing) > 0) {
        cat("        [WARN] Missing genes for combination:", paste(missing, collapse = ", "), "\n")
        return(NA)
    }

    tryCatch(
        {
            df <- data.frame(label = labels, mat[, gene_names, drop = FALSE], check.names = FALSE)

            # Build formula – backtick gene names to handle special characters (e.g. hyphens)
            gene_terms <- paste0("`", gene_names, "`")
            formula_str <- paste("label ~", paste(gene_terms, collapse = " + "))

            glm_fit <- glm(as.formula(formula_str), data = df, family = binomial(link = "logit"))
            preds <- predict(glm_fit, newdata = df, type = "response")

            roc_obj <- roc(labels, preds, quiet = TRUE, direction = "auto")
            as.numeric(auc(roc_obj))
        },
        error = function(e) {
            cat("        [ERROR] GLM failed:", conditionMessage(e), "\n")
            NA
        }
    )
}

################################################################################
# ── 4. Main Loop Over Conditions ─────────────────────────────────────────────
################################################################################

for (cond in conditions) {
    cat("\n################################################################\n")
    cat("  CONDITION:", cond$label, "\n")
    cat("################################################################\n\n")

    # Create condition-level output directory
    cond_out <- file.path(out_dir, cond$label)
    if (!dir.exists(cond_out)) dir.create(cond_out, recursive = TRUE)

    # ── 4a. Read top genes & combinations ────────────────────────────────────
    cat("  >> Reading gene statistics:", cond$gene_stats_file, "\n")
    gene_stats <- read.csv(cond$gene_stats_file, check.names = FALSE)
    top_genes <- gene_stats$Gene
    cat("     Top genes:", paste(top_genes, collapse = ", "), "\n\n")

    cat("  >> Reading GLM combinations:", cond$glm_combos_file, "\n")
    combos_df <- read.csv(cond$glm_combos_file, check.names = FALSE)
    combo_strings <- combos_df$Combination
    cat("     Number of combinations:", length(combo_strings), "\n\n")

    # ── 4b. Prepare the 3 validation datasets ───────────────────────────────
    cat("  >> Preparing GSE184316 dataset ...\n")
    ds_gse184316 <- prepare_gse184316(
        counts_file = cond$gse184316_counts,
        info_file   = cond$gse184316_info,
        pos_label   = cond$gse184316_pos
    )

    cat("  >> Preparing GSE150910 Biopsy dataset ...\n")
    ds_biopsy <- prepare_gse150910_subset(
        counts_file = cond$gse150910_counts,
        sample_file = cond$gse150910_biopsy
    )

    cat("  >> Preparing GSE150910 Explant dataset ...\n")
    ds_explant <- prepare_gse150910_subset(
        counts_file = cond$gse150910_counts,
        sample_file = cond$gse150910_explant
    )

    # Collect all 3 datasets into a named list for easy iteration
    datasets <- list(
        GSE184316         = ds_gse184316,
        GSE150910_Biopsy  = ds_biopsy,
        GSE150910_Explant = ds_explant
    )

    # ── 4c. Single-Gene AUC ─────────────────────────────────────────────────
    cat("\n  >> Computing single-gene AUC ...\n")

    single_results <- data.frame(Gene = top_genes, stringsAsFactors = FALSE)

    for (ds_name in names(datasets)) {
        ds <- datasets[[ds_name]]
        auc_vals <- sapply(top_genes, function(g) {
            compute_single_gene_auc(g, ds$mat, ds$labels)
        })
        single_results[[paste0("AUC_", ds_name)]] <- auc_vals
    }

    single_out_file <- file.path(cond_out, paste0(cond$label, "_single_gene_AUC.csv"))
    write.csv(single_results, single_out_file, row.names = FALSE)
    cat("     Saved:", single_out_file, "\n")

    # Print summary
    cat("\n     --- Single Gene AUC Summary ---\n")
    print(single_results)
    cat("\n")

    # ── 4d. Combination-Gene AUC (GLM) ─────────────────────────────────────
    cat("  >> Computing combination AUC (GLM) ...\n")

    combo_results <- data.frame(
        Combination = combo_strings,
        NumGenes = combos_df$NumGenes,
        stringsAsFactors = FALSE
    )

    for (ds_name in names(datasets)) {
        ds <- datasets[[ds_name]]

        auc_vals <- sapply(combo_strings, function(cs) {
            # Parse the gene names from the combination string
            gene_names <- trimws(unlist(strsplit(cs, ",")))
            compute_combo_auc(gene_names, ds$mat, ds$labels)
        })

        combo_results[[paste0("AUC_", ds_name)]] <- auc_vals
    }

    combo_out_file <- file.path(cond_out, paste0(cond$label, "_combination_AUC.csv"))
    write.csv(combo_results, combo_out_file, row.names = FALSE)
    cat("     Saved:", combo_out_file, "\n")

    # Print first few rows
    cat("\n     --- Combination AUC Summary (first 10 rows) ---\n")
    print(head(combo_results, 10))
    cat("\n")
}

################################################################################
# ── 5. logFC Heatmaps ────────────────────────────────────────────────────────
#   For each condition, read the gene_statistics.csv, extract the per-dataset
#   log2FC columns (GSE184316, GSE150910.biopsy, GSE150910.explant), and plot a
#   heatmap with genes as rows and datasets as columns.
#   Colour scale: blue (negative) → white (0) → red (positive).
################################################################################

cat("\n============================================================\n")
cat("  Section 5: logFC Heatmaps\n")
cat("============================================================\n\n")

suppressPackageStartupMessages({
    library(pheatmap)
    library(RColorBrewer)
})

# log2FC columns present in each gene_statistics.csv
logfc_columns <- c("GSE184316", "GSE150910.biopsy", "GSE150910.explant")

for (cond in conditions) {
    cat("  >> Generating logFC heatmap for:", cond$label, "\n")

    # ── Read gene statistics ─────────────────────────────────────────────────
    stats_file <- cond$gene_stats_file
    if (!file.exists(stats_file)) {
        cat("     [WARN] File not found, skipping:", stats_file, "\n")
        next
    }

    stats_df <- read.csv(stats_file, check.names = FALSE)

    if (nrow(stats_df) == 0) {
        cat("     [WARN] No genes found. Skipping.\n")
        next
    }

    cat("     Total genes:", nrow(stats_df), "\n")

    # ── Build numeric matrix ─────────────────────────────────────────────────
    mat <- as.matrix(stats_df[, logfc_columns])
    rownames(mat) <- stats_df$Gene

    # Prettier column labels for the heatmap
    colnames(mat) <- c("GSE184316", "GSE150910\n(Biopsy)", "GSE150910\n(Explant)")

    # ── Save log2FC matrix as CSV ─────────────────────────────────────────────
    csv_mat <- mat
    colnames(csv_mat) <- c("GSE184316", "GSE150910.biopsy", "GSE150910.explant")
    csv_file <- file.path(
        file.path(base_dir, "results", "boxplots", cond$label),
        paste0(cond$label, "_logFC_matrix.csv")
    )
    write.csv(csv_mat, csv_file)
    cat("     Saved log2FC matrix:", csv_file, "\n")

    # ── Symmetric colour breaks centred at 0 ─────────────────────────────────
    max_abs <- max(abs(mat))
    n_breaks <- 101
    breaks <- seq(-max_abs, max_abs, length.out = n_breaks)

    # Blue → White → Red
    colors <- colorRampPalette(c(
        "#2166AC", "#67A9CF", "#D1E5F0",
        "white",
        "#FDDBC7", "#EF8A62", "#B2182B"
    ))(n_breaks - 1)

    # ── Output path ──────────────────────────────────────────────────────────
    heatmap_dir <- file.path(base_dir, "results", "boxplots", cond$label)
    if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir, recursive = TRUE)
    heatmap_file <- file.path(heatmap_dir, paste0(cond$label, "_logFC_heatmap.png"))

    # ── Determine appropriate figure height (tight) ─────────────────────────
    # Compute height based on content: title + col labels + cells + legend
    n_genes <- nrow(mat)
    fig_height <- max(2.5, 0.30 * n_genes + 1.8)
    fig_width <- 5.5

    # ── Draw heatmap (column labels on top, no cell numbers) ─────────────────
    #    pheatmap places column labels at the bottom by default. We capture the
    #    grob, insert a new row above the matrix for column labels, move them
    #    there, and remove their original row to eliminate white space.
    p <- pheatmap(
        mat,
        cluster_rows      = TRUE,
        cluster_cols      = FALSE, # keep dataset order fixed
        scale             = "none",
        color             = colors,
        breaks            = breaks,
        border_color      = "grey80",
        cellwidth         = 60,
        cellheight        = 18,
        fontsize_row      = 9,
        fontsize_col      = 10,
        angle_col         = 0, # horizontal column labels
        main              = paste0(cond$label, " – log2 Fold Change by Dataset"),
        display_numbers   = FALSE,
        legend            = TRUE,
        silent            = TRUE # return grob without drawing
    )

    # ── Reposition column labels from bottom to just above the matrix ────────
    gt <- p$gtable

    col_idx <- which(gt$layout$name == "col_names")
    mat_idx <- which(gt$layout$name == "matrix")

    if (length(col_idx) == 1 && length(mat_idx) == 1) {
        mat_top_row <- gt$layout$t[mat_idx]

        # Insert a new row above the matrix for column labels
        gt <- gtable::gtable_add_rows(gt,
            heights = grid::unit(1.2, "cm"),
            pos = mat_top_row - 1
        )

        # After insertion, all indices at or below mat_top_row shift down by 1.
        # Re-find the indices after the structural change.
        col_idx <- which(gt$layout$name == "col_names")
        mat_idx <- which(gt$layout$name == "matrix")

        # Move col_names into the newly inserted row (which is now at mat_top_row)
        gt$layout$t[col_idx] <- gt$layout$t[mat_idx] - 1
        gt$layout$b[col_idx] <- gt$layout$t[mat_idx] - 1
    }

    # ── Trim excess padding from the gtable ──────────────────────────────────
    # Collapse tiny padding rows/cols that pheatmap inserts (< 3 pt)
    for (i in seq_along(gt$heights)) {
        h_pt <- tryCatch(
            as.numeric(grid::convertUnit(gt$heights[i], "pt")),
            error = function(e) NA
        )
        if (!is.na(h_pt) && h_pt >= 0 && h_pt < 3) {
            gt$heights[i] <- grid::unit(0, "pt")
        }
    }
    for (i in seq_along(gt$widths)) {
        w_pt <- tryCatch(
            as.numeric(grid::convertUnit(gt$widths[i], "pt")),
            error = function(e) NA
        )
        if (!is.na(w_pt) && w_pt >= 0 && w_pt < 3) {
            gt$widths[i] <- grid::unit(0, "pt")
        }
    }

    # ── Save with tight bounding box ─────────────────────────────────────────
    png(heatmap_file, width = fig_width, height = fig_height, units = "in", res = 300)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
        width  = grid::unit(0.95, "npc"),
        height = grid::unit(0.95, "npc")
    ))
    grid::grid.draw(gt)
    grid::popViewport()
    dev.off()

    cat("     Saved heatmap:", heatmap_file, "\n\n")
}

cat("\n============================================================\n")
cat("  All conditions processed. Results saved to:\n")
cat("  ", out_dir, "\n")
cat("  logFC heatmaps saved to: ./results/boxplots/{condition}/\n")
cat("============================================================\n")
