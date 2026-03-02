##############################################################################
## 07_boxplot.R
##
## PURPOSE:
##   Generate publication-quality comparison boxplots for the selected
##   predictor genes identified by the LASSO/Ridge pipeline (06_lasso_refactored.R).
##   A Mann-Whitney U test (Wilcoxon rank-sum) is used to assess
##   significance between the two groups in each pairwise condition.
##
## INPUT:
##   ./results/lasso/<condition>/lasso_selected_genes.txt  (one gene per line)
##   ./train_data/<condition>_train_data.csv               (samples × genes + condition)
##   ./test_data/<condition>_test_data.csv
##
## OUTPUT:
##   ./results/boxplots/<condition>/train/  — train-set boxplots (PNG + PDF)
##   ./results/boxplots/<condition>/test/   — test-set  boxplots (PNG + PDF)
##   ./results/boxplots/<condition>/<condition>_gene_statistics.csv
##       Columns: Gene, Mean_<neg>_Train, Mean_<pos>_Train, Mean_<neg>_Test,
##                Mean_<pos>_Test, pval_Train, pval_Test, padj_Train (BH),
##                padj_Test (BH), log2FC_Train, log2FC_Test,
##                + columns from {condition}_Common_genes_same_direction_expr.csv:
##                  ng, mean_abs_logfc, GSE184316, GSE150910.biopsy,
##                  GSE150910.explant, metaFIN_pval, metaFIN_padj
##
## CONDITIONS:
##   1) cHP_Ctrl  — HP  vs Control
##   2) cHP_IPF   — HP  vs IPF
##   3) IPF_Ctrl  — IPF vs Control
##############################################################################


# ──────────────────────────────────────────────────────────────────────────────
# 0. SETUP
# ──────────────────────────────────────────────────────────────────────────────

cat("============================================================\n")
cat("  Boxplot Pipeline — Pairwise Gene Expression Comparisons\n")
cat("============================================================\n\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

cat("[INFO] Libraries loaded: ggplot2, ggpubr\n\n")


# ──────────────────────────────────────────────────────────────────────────────
# 1. CONDITION DEFINITIONS
# ──────────────────────────────────────────────────────────────────────────────
# Each entry maps a condition tag to its human-readable label,
# positive / negative group names, and file paths.

conditions <- list(
  cHP_Ctrl = list(
    label      = "cHP vs Control",
    positive   = "HP",
    negative   = "Control",
    train_csv  = "./train_data/cHP_Ctrl_train_data.csv",
    test_csv   = "./test_data/cHP_Ctrl_test_data.csv",
    gene_file  = "./results/lasso/cHP_Ctrl/lasso_selected_genes.txt"
  ),
  cHP_IPF = list(
    label      = "cHP vs IPF",
    positive   = "HP",
    negative   = "IPF",
    train_csv  = "./train_data/cHP_IPF_train_data.csv",
    test_csv   = "./test_data/cHP_IPF_test_data.csv",
    gene_file  = "./results/lasso/cHP_IPF/lasso_selected_genes.txt"
  ),
  IPF_Ctrl = list(
    label      = "IPF vs Control",
    positive   = "IPF",
    negative   = "Control",
    train_csv  = "./train_data/IPF_Ctrl_train_data.csv",
    test_csv   = "./test_data/IPF_Ctrl_test_data.csv",
    gene_file  = "./results/lasso/IPF_Ctrl/lasso_selected_genes.txt"
  )
)


# ──────────────────────────────────────────────────────────────────────────────
# 2. HELPER: ensure_dir
# ──────────────────────────────────────────────────────────────────────────────

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("  [DIR]  Created: %s\n", path))
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# 3. CORE FUNCTION: plot_gene_boxplot
# ──────────────────────────────────────────────────────────────────────────────
#' Create a beautiful comparison boxplot for a single gene and save as
#' PNG (600 DPI) and PDF.
#'
#' @param gene         Character — gene (column) name to plot
#' @param dataset      Data frame containing gene columns + "condition" column
#' @param group_pos    Character — positive group label  (e.g. "HP")
#' @param group_neg    Character — negative group label  (e.g. "Control")
#' @param cond_label   Character — human-readable condition label for the title
#' @param dataset_name Character — "Train" or "Test" (used in title & filename)
#' @param out_dir      Character — directory in which to save the plots

plot_gene_boxplot <- function(gene, dataset, group_pos, group_neg,
                              cond_label, dataset_name, out_dir) {
  # --- Build a tidy data frame ------------------------------------------
  expr_df <- data.frame(
    expression = dataset[[gene]],
    group = dataset$condition,
    stringsAsFactors = FALSE
  )

  # Keep only the two groups of interest and set factor order
  expr_df <- expr_df[expr_df$group %in% c(group_pos, group_neg), ]
  expr_df$group <- factor(expr_df$group, levels = c(group_neg, group_pos))

  # --- Colour palette ---------------------------------------------------
  #   Soft, distinguishable colours that look great on both screen & print.
  fill_cols <- c("#4A90D9", "#E06666") # blue / coral
  border_cols <- c("#2C5F9E", "#B83C3C") # darker borders
  names(fill_cols) <- c(group_neg, group_pos)
  names(border_cols) <- c(group_neg, group_pos)

  # --- Build the plot ---------------------------------------------------
  p <- ggplot(expr_df, aes(x = group, y = expression, fill = group)) +

    # Boxplot: semi-transparent fill, coloured border, no outlier symbols
    geom_boxplot(
      aes(colour = group),
      width = 0.50,
      alpha = 0.35,
      outlier.shape = NA,
      lwd = 0.7,
      fatten = 1.5 # median line thickness multiplier
    ) +

    # Jittered individual data points
    geom_jitter(
      aes(colour = group),
      width = 0.18,
      size = 1.8,
      alpha = 0.75,
      shape = 16
    ) +

    # Colour / fill scales
    scale_fill_manual(values = fill_cols) +
    scale_colour_manual(values = border_cols) +

    # Mann-Whitney U test (Wilcoxon rank-sum)
    stat_compare_means(
      comparisons = list(c(group_neg, group_pos)),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.02,
      bracket.size = 0.5,
      size = 5,
      vjust = 0.3
    ) +

    # Labels
    labs(
      title    = paste0(gene, "  —  ", cond_label),
      subtitle = paste0(dataset_name, " set  |  Mann-Whitney U test"),
      y        = paste0(gene, " expression"),
      x        = NULL
    ) +

    # Theme: clean, modern, publication-ready
    theme_classic(base_size = 14) +
    theme(
      # Title styling
      plot.title = element_text(
        face = "bold", size = 15,
        hjust = 0.5, margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        size = 10, hjust = 0.5,
        colour = "grey40",
        margin = margin(b = 10)
      ),
      # Axis
      axis.title.y = element_text(
        face = "bold", size = 12,
        margin = margin(r = 8)
      ),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.line = element_line(linewidth = 0.5),

      # Legend
      legend.position = "none", # self-evident from x-axis

      # Margins
      plot.margin = margin(15, 15, 10, 15)
    )

  # --- Save as PNG & PDF ------------------------------------------------
  safe_gene <- gsub("[^a-zA-Z0-9_.-]", "_", gene) # sanitise for filenames
  base_name <- paste0(safe_gene, "_", dataset_name)

  ggsave(
    filename = file.path(out_dir, paste0(base_name, ".png")),
    plot = p,
    width = 5, height = 5.5, units = "in", dpi = 600
  )
  ggsave(
    filename = file.path(out_dir, paste0(base_name, ".pdf")),
    plot = p,
    width = 5, height = 5.5, units = "in"
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. MAIN LOOP — iterate over each condition
# ══════════════════════════════════════════════════════════════════════════════

for (cond_tag in names(conditions)) {
  cond <- conditions[[cond_tag]]

  cat("\n")
  cat("############################################################\n")
  cat(sprintf("##  CONDITION: %s  (%s)\n", cond$label, cond_tag))
  cat("############################################################\n\n")

  # ── 4A. Read the selected gene list ───────────────────────────────────
  if (!file.exists(cond$gene_file)) {
    cat(sprintf("  [SKIP] Gene file not found: %s\n", cond$gene_file))
    next
  }

  genes <- readLines(cond$gene_file)
  genes <- trimws(genes)
  genes <- genes[nchar(genes) > 0]
  cat(sprintf(
    "  [LOAD] %d selected genes from %s\n",
    length(genes), cond$gene_file
  ))
  cat(sprintf("         %s\n", paste(genes, collapse = ", ")))

  # ── 4B. Load training & test data ─────────────────────────────────────
  cat(sprintf("  [LOAD] Training data : %s\n", cond$train_csv))
  train_data <- read.csv(cond$train_csv, row.names = 1, check.names = FALSE)

  cat(sprintf("  [LOAD] Test data     : %s\n", cond$test_csv))
  test_data <- read.csv(cond$test_csv, row.names = 1, check.names = FALSE)

  cat(sprintf(
    "         Train samples: %d  |  Test samples: %d\n",
    nrow(train_data), nrow(test_data)
  ))

  # ── 4C. Verify genes exist in both datasets ───────────────────────────
  available <- intersect(genes, colnames(train_data))
  if (length(available) < length(genes)) {
    dropped <- setdiff(genes, available)
    cat(sprintf(
      "  [WARN] Gene(s) not found in data (dropped): %s\n",
      paste(dropped, collapse = ", ")
    ))
    genes <- available
  }

  if (length(genes) == 0) {
    cat("  [SKIP] No valid genes remaining. Skipping condition.\n")
    next
  }

  # ── 4D. Create output directories ─────────────────────────────────────
  out_train <- file.path(".", "results", "boxplots", cond_tag, "train")
  out_test <- file.path(".", "results", "boxplots", cond_tag, "test")
  ensure_dir(out_train)
  ensure_dir(out_test)

  # ── 4E. Generate boxplots for each gene ───────────────────────────────
  cat(sprintf(
    "\n  [PLOT] Generating boxplots for %d genes...\n",
    length(genes)
  ))

  for (gene in genes) {
    # ----- Training set -----
    plot_gene_boxplot(
      gene         = gene,
      dataset      = train_data,
      group_pos    = cond$positive,
      group_neg    = cond$negative,
      cond_label   = cond$label,
      dataset_name = "Train",
      out_dir      = out_train
    )

    # ----- Test set -----
    plot_gene_boxplot(
      gene         = gene,
      dataset      = test_data,
      group_pos    = cond$positive,
      group_neg    = cond$negative,
      cond_label   = cond$label,
      dataset_name = "Test",
      out_dir      = out_test
    )

    cat(sprintf("    ✓ %s\n", gene))
  }

  # ── 4F. Build per-gene statistics CSV ──────────────────────────────────
  cat("\n  ── Section 4F: Per-Gene Statistics CSV ────────────────────\n")

  # Pre-allocate vectors for each column
  mean_neg_train <- numeric(length(genes))
  mean_pos_train <- numeric(length(genes))
  mean_neg_test <- numeric(length(genes))
  mean_pos_test <- numeric(length(genes))
  pval_train <- numeric(length(genes))
  pval_test <- numeric(length(genes))

  for (i in seq_along(genes)) {
    gene <- genes[i]

    # ── Training set ──
    vals_neg_tr <- train_data[[gene]][train_data$condition == cond$negative]
    vals_pos_tr <- train_data[[gene]][train_data$condition == cond$positive]
    mean_neg_train[i] <- mean(vals_neg_tr, na.rm = TRUE)
    mean_pos_train[i] <- mean(vals_pos_tr, na.rm = TRUE)
    pval_train[i] <- wilcox.test(vals_pos_tr, vals_neg_tr)$p.value

    # ── Test set ──
    vals_neg_te <- test_data[[gene]][test_data$condition == cond$negative]
    vals_pos_te <- test_data[[gene]][test_data$condition == cond$positive]
    mean_neg_test[i] <- mean(vals_neg_te, na.rm = TRUE)
    mean_pos_test[i] <- mean(vals_pos_te, na.rm = TRUE)
    pval_test[i] <- wilcox.test(vals_pos_te, vals_neg_te)$p.value
  }

  # BH-adjusted p-values (across all genes within this condition)
  padj_train <- p.adjust(pval_train, method = "BH")
  padj_test <- p.adjust(pval_test, method = "BH")

  # log2 Fold Change: log2(positive_mean / negative_mean)
  log2fc_train <- log2(mean_pos_train / mean_neg_train)
  log2fc_test <- log2(mean_pos_test / mean_neg_test)

  # Assemble data frame
  stats_df <- data.frame(
    Gene = genes,
    Mean_Train = mean_neg_train,
    Mean_Train2 = mean_pos_train,
    Mean_Test = mean_neg_test,
    Mean_Test2 = mean_pos_test,
    pval_Train = pval_train,
    pval_Test = pval_test,
    padj_Train_BH = padj_train,
    padj_Test_BH = padj_test,
    log2FC_Train = log2fc_train,
    log2FC_Test = log2fc_test,
    stringsAsFactors = FALSE
  )

  # Rename mean columns to include the actual group labels
  colnames(stats_df)[colnames(stats_df) == "Mean_Train"] <-
    paste0("Mean_", cond$negative, "_Train")
  colnames(stats_df)[colnames(stats_df) == "Mean_Train2"] <-
    paste0("Mean_", cond$positive, "_Train")
  colnames(stats_df)[colnames(stats_df) == "Mean_Test"] <-
    paste0("Mean_", cond$negative, "_Test")
  colnames(stats_df)[colnames(stats_df) == "Mean_Test2"] <-
    paste0("Mean_", cond$positive, "_Test")

  # ── 4G. Merge columns from Common_genes_same_direction_expr CSV ──────
  common_expr_file <- file.path(
    ".", "data",
    paste0(cond_tag, "_Common_genes_same_direction_expr.csv")
  )

  if (file.exists(common_expr_file)) {
    cat(sprintf(
      "\n  [LOAD] Loading common-gene expression data: %s\n",
      common_expr_file
    ))
    common_expr <- read.csv(common_expr_file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    # Match stats_df genes to Gene_symbol in the common expression file
    matched_rows <- match(stats_df$Gene, common_expr$Gene_symbol)

    # Identify columns to add (everything except Gene_symbol itself)
    extra_cols <- setdiff(colnames(common_expr), "Gene_symbol")
    cat(sprintf(
      "  [INFO] Adding %d columns from common expression file: %s\n",
      length(extra_cols), paste(extra_cols, collapse = ", ")
    ))

    # Bind matched columns to stats_df
    for (col in extra_cols) {
      stats_df[[col]] <- common_expr[[col]][matched_rows]
    }

    n_matched <- sum(!is.na(matched_rows))
    n_missing <- sum(is.na(matched_rows))
    cat(sprintf(
      "  [INFO] Matched %d / %d genes (missing: %d)\n",
      n_matched, length(matched_rows), n_missing
    ))
    if (n_missing > 0) {
      missing_genes <- stats_df$Gene[is.na(matched_rows)]
      cat(sprintf(
        "  [WARN] Genes not found in common expression file: %s\n",
        paste(missing_genes, collapse = ", ")
      ))
    }
  } else {
    cat(sprintf(
      "\n  [WARN] Common expression file not found: %s — skipping merge.\n",
      common_expr_file
    ))
  }

  # Sort by adjusted p-value (train) ascending
  stats_df <- stats_df[order(stats_df$padj_Train_BH), ]

  # Print summary
  cat("\n  [TABLE] Gene statistics summary:\n")
  print(stats_df, row.names = FALSE, digits = 4)

  # Save CSV
  stats_csv <- file.path(
    ".", "results", "boxplots", cond_tag,
    paste0(cond_tag, "_gene_statistics.csv")
  )
  write.csv(stats_df, stats_csv, row.names = FALSE)
  cat(sprintf("\n  [SAVE] Statistics CSV → %s\n", stats_csv))

  cat(sprintf(
    "\n  [DONE] %d genes × 2 datasets → %s\n",
    length(genes),
    file.path(".", "results", "boxplots", cond_tag)
  ))
}


# ──────────────────────────────────────────────────────────────────────────────
# 5. FINISHED
# ──────────────────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat("  ALL BOXPLOTS GENERATED SUCCESSFULLY\n")
cat("  Results saved in ./results/boxplots/\n")
cat("============================================================\n")
