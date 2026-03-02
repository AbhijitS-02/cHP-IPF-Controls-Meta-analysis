##############################################################################
##  03b_subset_common_genes.R
##
##  PURPOSE:
##    Subset each batch-corrected expression matrix to retain ONLY the genes
##    that appear in the corresponding Common_genes_same_direction_expr CSV.
##    This ensures downstream analyses (mRMR, LASSO, boxplots, AUC-ROC) use
##    a consistent gene universe and avoids NAs when merging back log2FC
##    columns later.
##
##  INPUT:
##    ./results/Corrected_matrix/Merged_{tag}_BatchCorrected_Expression.csv
##    ./{condition}_Common_genes_same_direction_expr.csv
##
##  OUTPUT:
##    ./results/Subsetted_matrix/Merged_{tag}_Subsetted_Expression.csv
##    (original batch-corrected files are NOT modified)
##
##  CONDITIONS (mapping between filename tags):
##    cHP_Ctrl  →  cHPvsCtrl
##    cHP_IPF   →  cHPvsIPF
##    IPF_Ctrl  →  IPFvsCtrl
##
##############################################################################

cat("============================================================\n")
cat("  03b_subset_common_genes.R\n")
cat("  Subset batch-corrected matrices to common-direction genes\n")
cat("============================================================\n\n")

# ── Condition mapping ────────────────────────────────────────────────────────
# Each entry links:
#   condition  : the short tag used in Common_genes filenames
#   matrix_tag : the tag used in Merged_*_BatchCorrected filenames
conditions <- list(
    list(condition = "cHP_Ctrl", matrix_tag = "cHPvsCtrl"),
    list(condition = "cHP_IPF", matrix_tag = "cHPvsIPF"),
    list(condition = "IPF_Ctrl", matrix_tag = "IPFvsCtrl")
)

corrected_dir <- file.path(".", "results", "Corrected_matrix")
output_dir <- file.path(".", "results", "Subsetted_corrected_matrix")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ── Main loop ────────────────────────────────────────────────────────────────
for (cond in conditions) {
    cat(sprintf(
        "\n── Condition: %s ──────────────────────────────────\n",
        cond$condition
    ))

    # --- File paths -------------------------------------------------------
    common_file <- file.path(
        ".", paste0(cond$condition, "_Common_genes_same_direction_expr.csv")
    )
    matrix_file <- file.path(
        corrected_dir,
        paste0("Merged_", cond$matrix_tag, "_BatchCorrected_Expression.csv")
    )

    # --- Check existence --------------------------------------------------
    if (!file.exists(common_file)) {
        cat("  [SKIP] Common-genes file not found:", common_file, "\n")
        next
    }
    if (!file.exists(matrix_file)) {
        cat("  [SKIP] Batch-corrected matrix not found:", matrix_file, "\n")
        next
    }

    # --- Read common genes ------------------------------------------------
    common_df <- read.csv(common_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    common_genes <- common_df$Gene_symbol
    cat(sprintf(
        "  [LOAD] Common genes file : %s  (%d genes)\n",
        common_file, length(common_genes)
    ))

    # --- Read batch-corrected matrix --------------------------------------
    expr_mat <- read.csv(matrix_file, row.names = 1, check.names = FALSE)
    cat(sprintf(
        "  [LOAD] Expression matrix : %s  (%d genes × %d samples)\n",
        matrix_file, nrow(expr_mat), ncol(expr_mat)
    ))

    # --- Subset to common genes -------------------------------------------
    shared <- intersect(rownames(expr_mat), common_genes)
    missing_from_matrix <- setdiff(common_genes, rownames(expr_mat))

    cat(sprintf(
        "  [INFO] Shared genes : %d / %d common genes\n",
        length(shared), length(common_genes)
    ))

    if (length(missing_from_matrix) > 0) {
        cat(sprintf(
            "  [WARN] %d common genes NOT in expression matrix (skipped)\n",
            length(missing_from_matrix)
        ))
    }

    if (length(shared) == 0) {
        cat("  [ERROR] No shared genes — skipping. Check gene symbol format.\n")
        next
    }

    expr_sub <- expr_mat[shared, , drop = FALSE]
    cat(sprintf("  [SUBSET] %d → %d genes\n", nrow(expr_mat), nrow(expr_sub)))

    # --- Save to separate directory (originals untouched) -----------------
    out_file <- file.path(
        output_dir,
        paste0("Merged_", cond$matrix_tag, "_Subsetted_Expression.csv")
    )
    write.csv(expr_sub, out_file)
    cat(sprintf("  [SAVE] %s\n", out_file))
}

cat("\n============================================================\n")
cat("  Done. Subsetted matrices saved to:\n")
cat(sprintf("   %s\n", output_dir))
cat("  Original batch-corrected files are unchanged.\n")
cat("============================================================\n")
