##############################################################################
## 06_lasso_refactored.R
##
## PURPOSE:
##   Build and evaluate Ridge regression + GLM classification models for
##   three pairwise comparisons of gene-expression data:
##     1) cHP vs Control
##     2) cHP vs IPF
##     3) IPF vs Control
##
## WORKFLOW (per condition):
##   A. Load training / test CSV and mRMR gene list; filter to mRMR genes,
##      scale features (test scaled using train mean & SD)
##   B. Incremental gene sweep – Ridge regression (alpha = 0, 10-fold CV)
##      on top-k mRMR genes (k = 2..N) → AUC & Accuracy vs # genes plot
##   C. Interactive pause: user reviews sweep results and enters optimal_k
##      (the number of top mRMR genes to carry forward)
##   D. Per-gene ROC-AUC on training data (panel + individual plots)
##   E. Per-gene ROC-AUC on test data (panel + individual plots)
##   F. Combinatorial search (2–5 genes) using unpenalised logistic
##      regression (glm) – AUC & Accuracy
##      → Top-10 combination ROC plots for both train & test sets
##
## INPUT:
##   ./train_data/<condition>_train_data.csv   (samples × genes + condition)
##   ./test_data/<condition>_test_data.csv
##   ./train_data/<condition>_top_mrmr_genes.txt  (one gene per line)
##
## OUTPUT:
##   ./results/lasso/<condition>/  — all plots and CSV results
##############################################################################


# ──────────────────────────────────────────────────────────────────────────────
# 0. SETUP: Load required libraries
# ──────────────────────────────────────────────────────────────────────────────

cat("============================================================\n")
cat("  LASSO Modelling Pipeline – 3 Pairwise Comparisons\n")
cat("============================================================\n\n")

# Suppress package startup messages for cleaner console output
suppressPackageStartupMessages({
    library(tidyverse)
    library(glmnet)
    library(ggplot2)
    library(reshape2)
    library(pROC)
})

cat("[INFO] All libraries loaded successfully.\n\n")


# ──────────────────────────────────────────────────────────────────────────────
# 1. DEFINE CONDITIONS & FILE PATHS
# ──────────────────────────────────────────────────────────────────────────────
# Each element is a named list that maps the condition tag to:
#   - label       : A human-readable comparison name
#   - positive    : The class encoded as 1 (the "case" group)
#   - negative    : The class encoded as 0 (the "control / reference" group)
#   - train_csv   : Path to the training CSV file
#   - test_csv    : Path to the test CSV file
#   - mrmr_file   : Path to the mRMR gene list (.txt, one gene per line)

conditions <- list(
    cHP_Ctrl = list(
        label      = "cHP vs Control",
        positive   = "HP",
        negative   = "Control",
        train_csv  = "./train_data/cHP_Ctrl_train_data.csv",
        test_csv   = "./test_data/cHP_Ctrl_test_data.csv",
        mrmr_file  = "./train_data/cHP_Ctrl_top_mrmr_genes.txt"
    ),
    cHP_IPF = list(
        label      = "cHP vs IPF",
        positive   = "HP",
        negative   = "IPF",
        train_csv  = "./train_data/cHP_IPF_train_data.csv",
        test_csv   = "./test_data/cHP_IPF_test_data.csv",
        mrmr_file  = "./train_data/cHP_IPF_top_mrmr_genes.txt"
    ),
    IPF_Ctrl = list(
        label      = "IPF vs Control",
        positive   = "IPF",
        negative   = "Control",
        train_csv  = "./train_data/IPF_Ctrl_train_data.csv",
        test_csv   = "./test_data/IPF_Ctrl_test_data.csv",
        mrmr_file  = "./train_data/IPF_Ctrl_top_mrmr_genes.txt"
    )
)


# ──────────────────────────────────────────────────────────────────────────────
# 2. HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

#' ensure_dir
#' Create a directory (and parents) if it does not already exist.
ensure_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        cat(sprintf("  [DIR]  Created output directory: %s\n", path))
    }
}

#' evaluate_combination
#' Fit an unpenalized logistic regression (glm) on a given gene subset
#' and return AUC, accuracy, and the ROC object.
#'
#' Rationale: The input gene combinations have already been filtered for
#' mutual independence via the mRMR algorithm. Applying an L1/L2 penalty
#' at this stage would introduce artificial coefficient shrinkage and
#' prevent an accurate evaluation of each combination's joint predictive
#' power. Standard logistic regression is therefore the appropriate model.
#'
#' @param gene_combo  Character vector of gene names
#' @param X_train_sc  Scaled training matrix
#' @param y_train_b   Binary training labels (0/1)
#' @param X_test_sc   Scaled test matrix
#' @param y_test_b    Binary test labels (0/1)
#' @return Named list with auc, acc, roc, pred_train (training predictions), genes
evaluate_combination <- function(gene_combo, X_train_sc, y_train_b,
                                 X_test_sc, y_test_b) {
    # Build data frames for glm (response + predictors)
    train_df <- as.data.frame(X_train_sc[, gene_combo, drop = FALSE])
    train_df$y <- y_train_b

    test_df <- as.data.frame(X_test_sc[, gene_combo, drop = FALSE])

    # Fit unpenalized logistic regression
    glm_fit <- glm(y ~ ., data = train_df, family = binomial(link = "logit"))

    # ── Test-set predictions ──
    pred_prob_test <- predict(glm_fit, newdata = test_df, type = "response")
    roc_obj <- roc(y_test_b, as.vector(pred_prob_test), quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    pred_lab <- ifelse(pred_prob_test > 0.5, 1, 0)
    acc_val <- mean(pred_lab == y_test_b)

    # ── Training-set predictions (for train ROC later) ──
    pred_prob_train <- predict(glm_fit, newdata = train_df, type = "response")
    roc_obj_train <- roc(y_train_b, as.vector(pred_prob_train), quiet = TRUE)

    return(list(
        auc       = auc_val,
        acc       = acc_val,
        roc       = roc_obj,
        roc_train = roc_obj_train,
        genes     = gene_combo
    ))
}


#' plot_single_roc
#' Save a single ROC curve as a PNG with AUC + 95 % CI annotation.
plot_single_roc <- function(roc_obj, title, filepath, line_col = "red") {
    auc_val <- as.numeric(auc(roc_obj))
    ci_vals <- ci.auc(roc_obj)
    ci_lo <- round(ci_vals[1], 3)
    ci_hi <- round(ci_vals[3], 3)

    png(filename = filepath, width = 8, height = 8, units = "in", res = 300)
    plot(roc_obj,
        col = line_col, lwd = 3,
        xlab = "False Positive Rate",
        ylab = "True Positive Rate",
        main = title,
        cex.main = 1.0, cex.lab = 1.4, cex.axis = 1.2
    )
    text(0.6, 0.1,
        paste0(
            "AUC = ", round(auc_val, 3),
            "\n(95% CI: ", ci_lo, " - ", ci_hi, ")"
        ),
        cex = 1.2, font = 2, col = "black"
    )
    dev.off()
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. MAIN LOOP — iterate over each condition
# ══════════════════════════════════════════════════════════════════════════════

for (cond_tag in names(conditions)) {
    cond <- conditions[[cond_tag]]

    cat("\n")
    cat("############################################################\n")
    cat(sprintf("##  CONDITION: %s  (%s)\n", cond$label, cond_tag))
    cat("############################################################\n\n")

    # ------------------------------------------------------------------
    # 3A. Create output directory
    # ------------------------------------------------------------------
    out_dir <- file.path(".", "results", "lasso", cond_tag)
    ensure_dir(out_dir)

    # ------------------------------------------------------------------
    # 3B. Load data
    # ------------------------------------------------------------------
    cat(sprintf("  [LOAD] Training data : %s\n", cond$train_csv))
    train_raw <- read.csv(cond$train_csv, row.names = 1, check.names = FALSE)

    cat(sprintf("  [LOAD] Test data     : %s\n", cond$test_csv))
    test_raw <- read.csv(cond$test_csv, row.names = 1, check.names = FALSE)

    cat(sprintf("  [LOAD] mRMR genes    : %s\n", cond$mrmr_file))
    mRMR_genes <- readLines(cond$mrmr_file)
    mRMR_genes <- trimws(mRMR_genes) # remove any whitespace
    mRMR_genes <- mRMR_genes[nchar(mRMR_genes) > 0] # drop empty lines

    cat(sprintf("         → %d mRMR genes loaded\n", length(mRMR_genes)))
    cat(sprintf(
        "         → Training samples: %d | Test samples: %d\n",
        nrow(train_raw), nrow(test_raw)
    ))

    # ------------------------------------------------------------------
    # 3C. Prepare feature matrices & binary labels
    # ------------------------------------------------------------------
    # Keep only the mRMR genes that actually appear in the data columns
    available_genes <- intersect(mRMR_genes, colnames(train_raw))
    if (length(available_genes) < length(mRMR_genes)) {
        dropped <- setdiff(mRMR_genes, available_genes)
        cat(sprintf(
            "  [WARN] %d mRMR gene(s) not found in data and dropped: %s\n",
            length(dropped), paste(dropped, collapse = ", ")
        ))
        mRMR_genes <- available_genes
    }

    # Build feature data frames using only the mRMR genes
    X_train <- train_raw[, mRMR_genes, drop = FALSE]
    X_test <- test_raw[, mRMR_genes, drop = FALSE]

    # Binary labels: positive class → 1, negative class → 0
    y_train_fac <- factor(train_raw$condition,
        levels = c(cond$negative, cond$positive)
    )
    y_test_fac <- factor(test_raw$condition,
        levels = c(cond$negative, cond$positive)
    )

    y_train_bin <- ifelse(train_raw$condition == cond$positive, 1, 0)
    y_test_bin <- ifelse(test_raw$condition == cond$positive, 1, 0)

    cat(sprintf(
        "  [INFO] Binary encoding: '%s' = 1 (positive), '%s' = 0 (negative)\n",
        cond$positive, cond$negative
    ))
    cat(sprintf(
        "         Train class balance: %s=%d, %s=%d\n",
        cond$positive, sum(y_train_bin == 1),
        cond$negative, sum(y_train_bin == 0)
    ))
    cat(sprintf(
        "         Test  class balance: %s=%d, %s=%d\n",
        cond$positive, sum(y_test_bin == 1),
        cond$negative, sum(y_test_bin == 0)
    ))

    # ------------------------------------------------------------------
    # 3D. Scale features (test scaled using training parameters)
    # ------------------------------------------------------------------
    cat("  [SCALE] Scaling features (z-score; test scaled on train params)...\n")

    X_train_scaled <- scale(X_train)
    X_test_scaled <- scale(X_test,
        center = attr(X_train_scaled, "scaled:center"),
        scale  = attr(X_train_scaled, "scaled:scale")
    )

    # ================================================================
    # SECTION A: INCREMENTAL GENE SWEEP (top-k mRMR genes)
    # ================================================================
    cat("\n  ── Section A: Incremental Gene Sweep ──────────────────────\n")

    gene_counts <- 2:length(mRMR_genes)
    auc_list <- numeric(length(mRMR_genes))
    acc_list <- numeric(length(mRMR_genes))
    selected_gene_sets <- vector("list", length(mRMR_genes))

    for (k in gene_counts) {
        genes_k <- mRMR_genes[1:k]
        X_tr_sub <- X_train_scaled[, genes_k, drop = FALSE]
        X_te_sub <- X_test_scaled[, genes_k, drop = FALSE]

        # Ridge regression (alpha = 0) with 10-fold cross-validation
        # Ridge keeps all features in the model (no variable selection)
        # while shrinking coefficients to handle multicollinearity.
        set.seed(123)
        cv_fit <- cv.glmnet(
            as.matrix(X_tr_sub), y_train_bin,
            alpha = 0, family = "binomial",
            type.measure = "auc", nfolds = 10
        )

        # Predict on test set
        pred_prob <- predict(cv_fit,
            newx = as.matrix(X_te_sub),
            s = "lambda.min", type = "response"
        )

        # AUC
        roc_obj <- roc(y_test_bin, as.vector(pred_prob), quiet = TRUE)
        auc_list[k] <- as.numeric(auc(roc_obj))

        # Accuracy
        pred_label <- ifelse(pred_prob > 0.5, 1, 0)
        acc_list[k] <- mean(pred_label == y_test_bin)

        # Store the gene list for this k
        selected_gene_sets[[k]] <- genes_k

        if (k %% 50 == 0 || k == length(mRMR_genes)) {
            cat(sprintf(
                "    k = %3d genes | AUC = %.3f | Acc = %.3f\n",
                k, auc_list[k], acc_list[k]
            ))
        }
    }

    # Combine results into a tidy data frame
    performance_df <- data.frame(
        NumGenes = gene_counts,
        AUC      = auc_list[gene_counts],
        Accuracy = acc_list[gene_counts]
    )

    # Save performance table
    perf_csv <- file.path(out_dir, "gene_sweep_performance.csv")
    write.csv(performance_df, perf_csv, row.names = FALSE)
    cat(sprintf("  [SAVE] Gene sweep performance → %s\n", perf_csv))

    # ── Plot: AUC & Accuracy vs. Number of Genes ──
    sweep_plot <- file.path(out_dir, "LASSO_AUC_and_Accuracy_plot.png")
    png(filename = sweep_plot, height = 6, width = 8, units = "in", res = 300)

    plot(performance_df$NumGenes, performance_df$AUC,
        type = "b", pch = 19, col = "blue",
        xlab = "Number of Genes", ylab = "AUC", ylim = c(0, 1),
        main = paste0(cond$label, " — LASSO: AUC & Accuracy vs. # Genes")
    )
    grid()

    par(new = TRUE)
    plot(performance_df$NumGenes, performance_df$Accuracy,
        type = "b", pch = 17, col = "red",
        axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1)
    )
    axis(side = 4, col = "red", col.axis = "red")
    mtext("Accuracy", side = 4, line = 3, col = "red")

    legend("bottomright",
        legend = c("AUC", "Accuracy"),
        col = c("blue", "red"), pch = c(19, 17)
    )
    dev.off()
    cat(sprintf("  [PLOT] Gene sweep plot → %s\n", sweep_plot))

    # ================================================================
    # SECTION B: OPTIMAL GENE SET SELECTION  (interactive)
    # ================================================================
    # The pipeline pauses here so the user can inspect:
    #   1. The AUC & Accuracy vs. # Genes plot  (sweep_plot)
    #   2. The gene_sweep_performance.csv table  (perf_csv)
    # and then manually choose optimal_k – the number of top mRMR genes
    # to carry forward into the per-gene ROC and combination analyses.

    cat("\n  ── Section B: Optimal Gene Set Selection (interactive) ───\n")
    cat("\n")
    cat("  ┌──────────────────────────────────────────────────────────┐\n")
    cat("  │  PIPELINE PAUSED — please review the outputs below and  │\n")
    cat("  │  choose the optimal number of genes (optimal_k).        │\n")
    cat("  └──────────────────────────────────────────────────────────┘\n")
    cat(sprintf("  [REVIEW] Gene sweep plot : %s\n", sweep_plot))
    cat(sprintf("  [REVIEW] Performance CSV : %s\n", perf_csv))
    cat(sprintf(
        "  [INFO]   Valid range for optimal_k: 2 .. %d\n",
        length(mRMR_genes)
    ))
    cat("\n")

    # ── Interactive input loop with validation ──
    repeat {
        user_input <- readline(
            prompt = sprintf(
                "  >> Enter optimal_k for '%s' (2-%d): ",
                cond$label, length(mRMR_genes)
            )
        )

        # Convert to integer; NA if non-numeric
        optimal_k <- suppressWarnings(as.integer(user_input))

        if (is.na(optimal_k)) {
            cat("  [ERROR] Not a valid integer. Please try again.\n")
            next
        }
        if (optimal_k < 2 || optimal_k > length(mRMR_genes)) {
            cat(sprintf(
                "  [ERROR] Value must be between 2 and %d. Please try again.\n",
                length(mRMR_genes)
            ))
            next
        }

        # Valid input — confirm and break
        cat(sprintf("  [OK]    optimal_k = %d accepted.\n", optimal_k))
        break
    }

    # Extract the top optimal_k mRMR genes
    lasso_genes <- selected_gene_sets[[optimal_k]]

    cat(sprintf(
        "  [INFO] Selected %d genes for downstream analysis:\n",
        optimal_k
    ))
    cat(sprintf("         %s\n", paste(lasso_genes, collapse = ", ")))

    # Save the selected gene list to file
    gene_file <- file.path(out_dir, "lasso_selected_genes.txt")
    writeLines(lasso_genes, gene_file)
    cat(sprintf("  [SAVE] Selected genes → %s\n", gene_file))

    # ================================================================
    # SECTION C: PER-GENE ROC-AUC ON TRAINING DATA
    # ================================================================
    cat("\n  ── Section C: Per-Gene ROC-AUC (Train) ───────────────────\n")

    train_roc_dir <- file.path(out_dir, "train_roc_individual")
    ensure_dir(train_roc_dir)

    # Accumulate per-gene AUC results
    auc_results_train <- data.frame(
        Gene = character(), AUC = numeric(),
        CI_lower = numeric(), CI_upper = numeric(),
        stringsAsFactors = FALSE
    )

    for (gene in lasso_genes) {
        roc_obj <- roc(y_train_bin, X_train[[gene]], quiet = TRUE)
        auc_val <- as.numeric(auc(roc_obj))
        ci_vals <- ci.auc(roc_obj)
        auc_results_train <- rbind(auc_results_train, data.frame(
            Gene     = gene,
            AUC      = auc_val,
            CI_lower = ci_vals[1],
            CI_upper = ci_vals[3]
        ))
    }

    # Sort descending by AUC
    auc_results_train <- auc_results_train[order(-auc_results_train$AUC), ]
    cat("  [TABLE] Per-gene AUC on training data:\n")
    print(auc_results_train, row.names = FALSE)

    # Save table
    write.csv(auc_results_train,
        file.path(out_dir, "per_gene_AUC_train.csv"),
        row.names = FALSE
    )

    # ── Panel ROC plot (all selected genes) ──
    n_genes <- length(lasso_genes)
    n_cols <- min(4, n_genes)
    n_rows <- ceiling(n_genes / n_cols)

    panel_file <- file.path(out_dir, "Predictor_genes_roc_train_panel.png")
    png(
        filename = panel_file, width = 2.5 * n_cols, height = 2.5 * n_rows,
        units = "in", res = 300
    )
    par(mfrow = c(n_rows, n_cols))

    for (gene in auc_results_train$Gene) {
        roc_obj <- roc(y_train_fac, X_train[[gene]], quiet = TRUE)
        plot(roc_obj,
            main = gene, col = "red",
            xlab = "FPR", ylab = "TPR"
        )
        legend("bottomright",
            legend = paste("AUC =", round(auc(roc_obj), 3)),
            cex = 0.7, bty = "n"
        )
    }
    dev.off()
    cat(sprintf("  [PLOT] Train ROC panel → %s\n", panel_file))

    # ── Individual ROC plots for each gene (train) ──
    for (gene in auc_results_train$Gene) {
        roc_obj <- roc(y_train_fac, X_train[[gene]], quiet = TRUE)
        fp <- file.path(train_roc_dir, paste0("Train_ROC_", gene, ".png"))
        plot_single_roc(roc_obj,
            title = paste0(gene, " (Train)"),
            filepath = fp, line_col = "red"
        )
    }
    cat(sprintf(
        "  [PLOT] %d individual train ROC plots → %s/\n",
        n_genes, train_roc_dir
    ))

    # ================================================================
    # SECTION D: PER-GENE ROC-AUC ON TEST DATA
    # ================================================================
    cat("\n  ── Section D: Per-Gene ROC-AUC (Test) ────────────────────\n")

    test_roc_dir <- file.path(out_dir, "test_roc_individual")
    ensure_dir(test_roc_dir)

    auc_results_test <- data.frame(
        Gene = character(), AUC = numeric(),
        CI_lower = numeric(), CI_upper = numeric(),
        stringsAsFactors = FALSE
    )

    for (gene in lasso_genes) {
        roc_obj <- roc(y_test_bin, X_test[[gene]], quiet = TRUE)
        auc_val <- as.numeric(auc(roc_obj))
        ci_vals <- ci.auc(roc_obj)
        auc_results_test <- rbind(auc_results_test, data.frame(
            Gene     = gene,
            AUC      = auc_val,
            CI_lower = ci_vals[1],
            CI_upper = ci_vals[3]
        ))
    }

    auc_results_test <- auc_results_test[order(-auc_results_test$AUC), ]
    cat("  [TABLE] Per-gene AUC on test data:\n")
    print(auc_results_test, row.names = FALSE)

    write.csv(auc_results_test,
        file.path(out_dir, "per_gene_AUC_test.csv"),
        row.names = FALSE
    )

    # ── Panel ROC plot (test) ──
    panel_file_test <- file.path(out_dir, "Predictor_genes_roc_test_panel.png")
    png(
        filename = panel_file_test, width = 2.5 * n_cols,
        height = 2.5 * n_rows, units = "in", res = 300
    )
    par(mfrow = c(n_rows, n_cols))

    for (gene in auc_results_test$Gene) {
        roc_obj <- roc(y_test_fac, X_test[[gene]], quiet = TRUE)
        plot(roc_obj,
            main = gene, col = "red",
            xlab = "FPR", ylab = "TPR"
        )
        legend("bottomright",
            legend = paste("AUC =", round(auc(roc_obj), 3)),
            cex = 0.7, bty = "n"
        )
    }
    dev.off()
    cat(sprintf("  [PLOT] Test ROC panel → %s\n", panel_file_test))

    # ── Individual ROC plots for each gene (test) ──
    for (gene in auc_results_test$Gene) {
        roc_obj <- roc(y_test_fac, X_test[[gene]], quiet = TRUE)
        fp <- file.path(test_roc_dir, paste0("Test_ROC_", gene, ".png"))
        plot_single_roc(roc_obj,
            title = paste0(gene, " (Test)"),
            filepath = fp, line_col = "red"
        )
    }
    cat(sprintf(
        "  [PLOT] %d individual test ROC plots → %s/\n",
        n_genes, test_roc_dir
    ))

    # ================================================================
    # SECTION E: GENE COMBINATIONS (2–5 genes from the selected set)
    # ================================================================
    # Uses unpenalized logistic regression (glm) instead of LASSO.
    # The mRMR-filtered genes are already selected for low redundancy,
    # so no further regularisation is needed — this gives a true
    # measure of each combination's joint predictive power.

    cat("\n  ── Section E: Combinatorial Gene Search (GLM) ────────────\n")

    combo_results <- data.frame(
        Combination = character(), NumGenes = integer(),
        AUC = numeric(), Accuracy = numeric(),
        stringsAsFactors = FALSE
    )
    top_combos_list <- list() # store test ROC objects
    top_combos_train_list <- list() # store train ROC objects

    # Limit the max combination size to avoid combinatorial explosion
    # if the selected gene set is very large
    max_combo_size <- min(5, length(lasso_genes))

    for (k in 2:max_combo_size) {
        combos <- combn(lasso_genes, k, simplify = FALSE)
        cat(sprintf(
            "    k=%d : evaluating %d combinations (glm)...\n",
            k, length(combos)
        ))

        for (combo in combos) {
            combo_name <- paste(combo, collapse = ", ")

            res <- evaluate_combination(
                combo, X_train_scaled, y_train_bin,
                X_test_scaled, y_test_bin
            )

            combo_results <- rbind(combo_results, data.frame(
                Combination = combo_name,
                NumGenes = k,
                AUC = res$auc,
                Accuracy = res$acc,
                stringsAsFactors = FALSE
            ))

            top_combos_list[[combo_name]] <- res$roc
            top_combos_train_list[[combo_name]] <- res$roc_train
        }
    }

    # Sort by AUC descending
    combo_results <- combo_results[order(-combo_results$AUC), ]

    cat("  [TABLE] Top 10 combinations by AUC:\n")
    print(head(combo_results, 10), row.names = FALSE)

    # Save full results
    combo_csv <- file.path(out_dir, "GLM_combinations_AUC_results.csv")
    write.csv(combo_results, combo_csv, row.names = FALSE)
    cat(sprintf("  [SAVE] Combination results → %s\n", combo_csv))

    # ── Top-10 ROC plots (Test) ──
    cat("\n  ── Section E.1: Top-10 Combination ROC Plots (Test) ─────\n")

    combo_test_dir <- file.path(out_dir, "Combination_Top10_Test_ROC")
    ensure_dir(combo_test_dir)

    top_10 <- head(combo_results, 10)

    for (i in seq_len(nrow(top_10))) {
        combo_name <- top_10$Combination[i]
        roc_obj <- top_combos_list[[combo_name]]
        safe_name <- gsub("[^a-zA-Z0-9_]", "_", combo_name)
        fp <- file.path(
            combo_test_dir,
            paste0("ROC_", i, "_", safe_name, ".png")
        )

        plot_single_roc(roc_obj,
            title = paste0("Combo ", i, " (Test):\n", combo_name),
            filepath = fp, line_col = "red"
        )
    }
    cat(sprintf("  [PLOT] 10 test ROC plots → %s/\n", combo_test_dir))

    # ── Top-10 ROC plots (Train) ──
    # Train ROC objects were already computed inside evaluate_combination()
    # and stored in top_combos_train_list, so no re-fitting is needed.
    cat("\n  ── Section E.2: Top-10 Combination ROC Plots (Train) ────\n")

    combo_train_dir <- file.path(out_dir, "Combination_Top10_Train_ROC")
    ensure_dir(combo_train_dir)

    for (i in seq_len(nrow(top_10))) {
        combo_name <- top_10$Combination[i]
        roc_obj_train <- top_combos_train_list[[combo_name]]

        safe_name <- gsub("[^a-zA-Z0-9_]", "_", combo_name)
        fp <- file.path(
            combo_train_dir,
            paste0("ROC_Train_", i, "_", safe_name, ".png")
        )

        plot_single_roc(roc_obj_train,
            title = paste0("Combo ", i, " (Train):\n", combo_name),
            filepath = fp, line_col = "blue"
        )
    }
    cat(sprintf("  [PLOT] 10 train ROC plots → %s/\n", combo_train_dir))

    cat(sprintf("\n  ✓ Condition '%s' complete.\n", cond$label))
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. DONE
# ──────────────────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat("  ALL CONDITIONS PROCESSED SUCCESSFULLY\n")
cat("  Results saved in ./results/lasso/\n")
cat("============================================================\n")
