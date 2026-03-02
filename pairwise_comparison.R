# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Merge metadata

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Define your input files (Adjust paths if necessary)
# We map the Score File to its corresponding Metadata File
comparisons <- list(
  list(
    name = "HP_vs_IPF",
    scores = "./results/xCell_results/Merged_scores_HP_IPF.csv",
    # REPLACE with your actual metadata filename for this comparison
    metadata = "./results/xCell_results/Merged_metadata_HP_IPF.csv",
    cond_col = "condition", # Column name in metadata for Group (e.g., 'condition')
    sample_col = "sample" # Column name in metadata for Sample ID (e.g., 'sample')
  ),
  list(
    name = "IPF_vs_Ctrl",
    scores = "./results/xCell_results/Merged_scores_IPF_Ctrl.csv",
    metadata = "./results/xCell_results/Merged_metadata_IPF_Ctrl.csv",
    cond_col = "condition",
    sample_col = "sample"
  ),
  list(
    name = "HP_vs_Ctrl",
    scores = "./results/xCell_results/Merged_scores_HP_Ctrl.csv",
    metadata = "./results/xCell_results/Merged_metadata_HP_Ctrl.csv",
    cond_col = "condition",
    sample_col = "sample"
  )
)

# Output directory
out_dir <- "./results/xCell_Stats_Plots/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================

run_xcell_analysis <- function(comp) {
  message(paste0("--- Processing: ", comp$name, " ---"))

  # 1. Load Data
  if (!file.exists(comp$scores)) stop(paste("File not found:", comp$scores))
  if (!file.exists(comp$metadata)) stop(paste("File not found:", comp$metadata))

  scores <- read.csv(comp$scores, row.names = 1, check.names = FALSE)
  meta <- read.csv(comp$metadata)

  # 2. Align Metadata to Scores
  # Ensure metadata rows match score columns
  # Filter metadata to only include samples present in the score file
  common_samples <- intersect(colnames(scores), meta[[comp$sample_col]])

  if (length(common_samples) == 0) stop("No matching sample IDs found between Scores and Metadata!")

  scores <- scores[, common_samples]
  meta <- meta[match(common_samples, meta[[comp$sample_col]]), ]

  # Check alignment
  if (!all(colnames(scores) == meta[[comp$sample_col]])) stop("Alignment failed!")

  # Get groups
  groups <- factor(meta[[comp$cond_col]])
  group_levels <- levels(groups)

  # 3. Initialize Results Dataframe
  results_df <- data.frame(
    Cell_Type = rownames(scores),
    Mean_Group1 = NA,
    Mean_Group2 = NA,
    P_Value = NA,
    FoldChange = NA,
    stringsAsFactors = FALSE
  )

  # Update column names to be specific (e.g., Mean_HP, Mean_Control)
  colnames(results_df)[2:3] <- paste0("Mean_", group_levels)

  # 4. Loop through Cell Types for Stats
  for (i in 1:nrow(scores)) {
    cell_type <- rownames(scores)[i]
    values <- as.numeric(scores[i, ])

    # Run Mann-Whitney U test
    # tryCatch handles cases where variance is 0 (all scores are identical)
    test_res <- tryCatch(
      {
        wilcox.test(values ~ groups, exact = FALSE)
      },
      error = function(e) {
        return(NULL)
      }
    )

    if (!is.null(test_res)) {
      results_df$P_Value[i] <- test_res$p.value
      # Store means
      means <- tapply(values, groups, mean)
      results_df[i, 2] <- means[1]
      results_df[i, 3] <- means[2]
      # Calculate Fold Change (Mean_Group2 / Mean_Group1)
      if (!is.na(means[1]) && means[1] != 0) {
        results_df$FoldChange[i] <- means[2] / means[1]
      }
    }
  }

  # 5. Apply FDR Correction (Benjamini-Hochberg)
  results_df$FDR <- p.adjust(results_df$P_Value, method = "BH")

  # Sort by P-value
  results_df <- results_df %>% arrange(P_Value)

  # Save Stats to CSV
  write.csv(results_df, paste0(out_dir, "Stats_", comp$name, ".csv"), row.names = FALSE)

  # ==============================================================================
  # 6. Generate PDF with Box Plots
  # ==============================================================================
  pdf_file <- paste0(out_dir, "BoxPlots_", comp$name, ".pdf")
  pdf(pdf_file, width = 6, height = 5)

  # Loop through ALL cell types (sorted by significance)
  for (i in 1:nrow(results_df)) {
    cell <- results_df$Cell_Type[i]
    pval <- signif(results_df$P_Value[i], 3)
    fdr <- signif(results_df$FDR[i], 3)

    # Prepare data for plotting
    plot_data <- data.frame(
      Score = as.numeric(scores[cell, ]),
      Condition = groups
    )

    # Define color based on significance
    title_color <- if (!is.na(fdr) && fdr < 0.05) "red" else "black"

    # Create Box Plot
    p <- ggplot(plot_data, aes(x = Condition, y = Score, fill = Condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) + # Boxplot
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) + # Jitter points
      theme_bw() +
      labs(
        title = cell,
        subtitle = paste0("p-value: ", pval, " | FDR: ", fdr),
        y = "xCell Enrichment Score",
        x = ""
      ) +
      theme(
        plot.title = element_text(face = "bold", color = title_color),
        legend.position = "none" # Hide legend as X-axis is labeled
      )

    print(p)
  }

  dev.off()
  message(paste("Saved PDF:", pdf_file))
}


# ==============================================================================
# BUBBLE/DOT HEATMAP - Combined across all comparisons
# ==============================================================================

generate_bubble_heatmap <- function(comparisons, out_dir) {
  message("--- Generating combined bubble heatmap ---")

  # Read and combine all stats files
  all_stats <- list()

  for (comp in comparisons) {
    stats_file <- paste0(out_dir, "Stats_", comp$name, ".csv")
    if (!file.exists(stats_file)) {
      message(paste("Stats file not found for", comp$name, "- skipping."))
      next
    }
    stats <- read.csv(stats_file)
    stats$Comparison <- comp$name
    all_stats[[comp$name]] <- stats
  }

  if (length(all_stats) == 0) {
    message("No stats files found. Cannot generate bubble heatmap.")
    return(NULL)
  }

  combined <- bind_rows(all_stats)

  # Replace NA FoldChange with 1 (no change) for plotting
  combined$FoldChange[is.na(combined$FoldChange)] <- 1

  # Create discrete FDR significance categories
  combined$FDR_Category <- cut(
    combined$FDR,
    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
    labels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "NS"),
    right = TRUE
  )
  # Handle NA FDR values as non-significant
  combined$FDR_Category[is.na(combined$FDR_Category)] <- "NS"

  # --- Remove cell types that are NS across ALL comparisons ---
  sig_cells <- combined %>%
    filter(FDR_Category != "NS") %>%
    pull(Cell_Type) %>%
    unique()
  combined <- combined %>% filter(Cell_Type %in% sig_cells)

  # --- Remove NS rows (no bubble shown for non-significant) ---
  combined <- combined %>% filter(FDR_Category != "NS")

  # Drop unused "NS" level from the factor
  combined$FDR_Category <- droplevels(combined$FDR_Category)

  # --- Rank cell types by number of comparisons they are significant in ---
  # Count how many comparisons each cell type appears in (after NS removal)
  sig_count <- combined %>%
    group_by(Cell_Type) %>%
    summarise(n_sig = n_distinct(Comparison), .groups = "drop") %>%
    arrange(desc(n_sig), Cell_Type) # Descending count, then alphabetical

  # Set Cell_Type as an ordered factor: sig in 3 → 2 → 1
  combined$Cell_Type <- factor(combined$Cell_Type, levels = sig_count$Cell_Type)

  # --- Compute log2 Fold Change for color ---
  combined$log2FC <- log2(combined$FoldChange)
  combined$log2FC[is.na(combined$log2FC) | is.infinite(combined$log2FC)] <- 0

  # Set comparison order
  combined$Comparison <- factor(combined$Comparison,
    levels = c("HP_vs_Ctrl", "HP_vs_IPF", "IPF_vs_Ctrl")
  )

  # Define discrete bubble sizes for FDR categories (larger = more significant)
  fdr_sizes <- c(
    "< 0.0001" = 8,
    "< 0.001"  = 6,
    "< 0.01"   = 4,
    "< 0.05"   = 2
  )

  # Determine symmetric color limits for log2FC
  max_abs_fc <- max(abs(combined$log2FC), na.rm = TRUE)
  if (max_abs_fc == 0) max_abs_fc <- 1

  # Build the publication-ready bubble plot (landscape: cell types on X-axis)
  p <- ggplot(combined, aes(x = Cell_Type, y = Comparison)) +
    geom_point(aes(size = FDR_Category, fill = log2FC),
      shape = 21, color = "black", stroke = 0.4
    ) +
    scale_fill_gradientn(
      colours = c(
        "#2166AC", "#4393C3", "#92C5DE",
        "white",
        "#F4A582", "#D6604D", "#B2182B"
      ),
      values = scales::rescale(
        c(
          -max_abs_fc, -max_abs_fc * 0.3, -max_abs_fc * 0.08,
          0,
          max_abs_fc * 0.08, max_abs_fc * 0.3, max_abs_fc
        ),
        to = c(0, 1)
      ),
      limits = c(-max_abs_fc, max_abs_fc),
      name = expression(log[2](FC))
    ) +
    scale_size_manual(
      values = fdr_sizes,
      name = "FDR",
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      # Panel — no grid
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      # Axes
      axis.text.x = element_text(
        size = 7.5, angle = 90, hjust = 1, vjust = 0.5,
        color = "black", face = "plain"
      ),
      axis.text.y = element_text(size = 10, color = "black", face = "bold"),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      # Legend
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.margin = margin(t = 5, b = 5),
      # Plot margins
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    ) +
    guides(
      fill = guide_colorbar(
        order = 1, barwidth = 8, barheight = 1,
        title.position = "top"
      ),
      size = guide_legend(
        order = 2, title.position = "top",
        override.aes = list(fill = "grey50")
      )
    )


  # Determine plot width dynamically based on number of cell types
  n_cells <- length(unique(combined$Cell_Type))
  plot_width <- max(10, n_cells * 0.35 + 4)
  plot_height <- 7.5

  # Save as landscape PDF
  pdf_file <- paste0(out_dir, "BubbleHeatmap_AllComparisons.pdf")
  ggsave(pdf_file,
    plot = p, width = plot_width, height = plot_height, limitsize = FALSE,
    dpi = 1200
  )
  message(paste("Saved bubble heatmap:", pdf_file))

  # Save as landscape PNG
  png_file <- paste0(out_dir, "BubbleHeatmap_AllComparisons.png")
  ggsave(png_file,
    plot = p, width = plot_width, height = plot_height, limitsize = FALSE,
    dpi = 1200
  )
  message(paste("Saved bubble heatmap:", png_file))

  return(p)
}


# ==============================================================================
# EXECUTE
# ==============================================================================

# Step 1: Run pairwise analyses for all comparisons
for (comp in comparisons) {
  if (file.exists(comp$scores) && file.exists(comp$metadata)) {
    run_xcell_analysis(comp)
  } else {
    message(paste("Skipping", comp$name, "- Metadata or Score file missing."))
  }
}

# Step 2: Generate a single combined bubble heatmap
generate_bubble_heatmap(comparisons, out_dir)


# ==============================================================================
# BATCH-CORRECTED ANALYSIS
# ==============================================================================
# Repeat the entire pairwise comparison pipeline using ComBat batch-corrected
# enrichment scores generated by 03_xCell.R.
# The corrected score files have cell types as rows and samples as columns,
# same structure as the raw merged scores but with platform effects removed.
# We reuse the SAME metadata files since sample identities are unchanged.
# NOTE: Corrected files only contain cell types that passed the zero-variance
# filter in 03_xCell.R, so the number of cell types may be slightly fewer.
# ==============================================================================

message("\n========== BATCH-CORRECTED ENRICHMENT SCORES ==========\n")

comparisons_corrected <- list(
  list(
    name = "HP_vs_IPF",
    scores = "./results/xCell_results/Corrected_scores/Corrected_scores_HP_IPF.csv",
    metadata = "./results/xCell_results/Merged_metadata_HP_IPF.csv",
    cond_col = "condition",
    sample_col = "sample"
  ),
  list(
    name = "IPF_vs_Ctrl",
    scores = "./results/xCell_results/Corrected_scores/Corrected_scores_IPF_Ctrl.csv",
    metadata = "./results/xCell_results/Merged_metadata_IPF_Ctrl.csv",
    cond_col = "condition",
    sample_col = "sample"
  ),
  list(
    name = "HP_vs_Ctrl",
    scores = "./results/xCell_results/Corrected_scores/Corrected_scores_HP_Ctrl.csv",
    metadata = "./results/xCell_results/Merged_metadata_HP_Ctrl.csv",
    cond_col = "condition",
    sample_col = "sample"
  )
)

# Output directory for corrected results
out_dir_corrected <- "./results/xCell_Stats_Plots_Corrected/"
dir.create(out_dir_corrected, showWarnings = FALSE, recursive = TRUE)

# Step 1: Run pairwise analyses for all corrected comparisons
for (comp in comparisons_corrected) {
  if (file.exists(comp$scores) && file.exists(comp$metadata)) {
    # Temporarily override out_dir used by run_xcell_analysis
    # We pass a modified comp with the corrected out_dir
    run_xcell_analysis_corrected <- function(comp, out_dir_corr) {
      message(paste0("--- Processing (Corrected): ", comp$name, " ---"))

      # 1. Load Data
      if (!file.exists(comp$scores)) stop(paste("File not found:", comp$scores))
      if (!file.exists(comp$metadata)) stop(paste("File not found:", comp$metadata))

      scores <- read.csv(comp$scores, row.names = 1, check.names = FALSE)
      meta <- read.csv(comp$metadata)

      # 2. Align Metadata to Scores
      common_samples <- intersect(colnames(scores), meta[[comp$sample_col]])

      if (length(common_samples) == 0) stop("No matching sample IDs found between Scores and Metadata!")

      scores <- scores[, common_samples]
      meta <- meta[match(common_samples, meta[[comp$sample_col]]), ]

      # Check alignment
      if (!all(colnames(scores) == meta[[comp$sample_col]])) stop("Alignment failed!")

      # Get groups
      groups <- factor(meta[[comp$cond_col]])
      group_levels <- levels(groups)

      # 3. Initialize Results Dataframe
      results_df <- data.frame(
        Cell_Type = rownames(scores),
        Mean_Group1 = NA,
        Mean_Group2 = NA,
        P_Value = NA,
        FoldChange = NA,
        stringsAsFactors = FALSE
      )

      colnames(results_df)[2:3] <- paste0("Mean_", group_levels)

      # 4. Loop through Cell Types for Stats
      for (i in 1:nrow(scores)) {
        cell_type <- rownames(scores)[i]
        values <- as.numeric(scores[i, ])

        test_res <- tryCatch(
          {
            wilcox.test(values ~ groups, exact = FALSE)
          },
          error = function(e) {
            return(NULL)
          }
        )

        if (!is.null(test_res)) {
          results_df$P_Value[i] <- test_res$p.value
          means <- tapply(values, groups, mean)
          results_df[i, 2] <- means[1]
          results_df[i, 3] <- means[2]
          if (!is.na(means[1]) && means[1] != 0) {
            results_df$FoldChange[i] <- means[2] / means[1]
          }
        }
      }

      # 5. Apply FDR Correction (Benjamini-Hochberg)
      results_df$FDR <- p.adjust(results_df$P_Value, method = "BH")

      # Sort by P-value
      results_df <- results_df %>% arrange(P_Value)

      # Save Stats to CSV
      write.csv(results_df, paste0(out_dir_corr, "Stats_", comp$name, ".csv"), row.names = FALSE)

      # 6. Generate PDF with Box Plots
      pdf_file <- paste0(out_dir_corr, "BoxPlots_", comp$name, ".pdf")
      pdf(pdf_file, width = 6, height = 5)

      for (i in 1:nrow(results_df)) {
        cell <- results_df$Cell_Type[i]
        pval <- signif(results_df$P_Value[i], 3)
        fdr <- signif(results_df$FDR[i], 3)

        plot_data <- data.frame(
          Score = as.numeric(scores[cell, ]),
          Condition = groups
        )

        title_color <- if (!is.na(fdr) && fdr < 0.05) "red" else "black"

        p <- ggplot(plot_data, aes(x = Condition, y = Score, fill = Condition)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.6) +
          geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
          theme_bw() +
          labs(
            title = cell,
            subtitle = paste0("p-value: ", pval, " | FDR: ", fdr),
            y = "xCell Enrichment Score (Corrected)",
            x = ""
          ) +
          theme(
            plot.title = element_text(face = "bold", color = title_color),
            legend.position = "none"
          )

        print(p)
      }

      dev.off()
      message(paste("Saved PDF:", pdf_file))
    }

    run_xcell_analysis_corrected(comp, out_dir_corrected)
  } else {
    message(paste("Skipping (Corrected)", comp$name, "- Score or Metadata file missing."))
  }
}

# Step 2: Generate a single combined bubble heatmap for corrected scores
generate_bubble_heatmap(comparisons_corrected, out_dir_corrected)
