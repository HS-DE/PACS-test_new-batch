long_metrics <- melt(
  as.data.table(strategy_summary),
  id.vars = c("strategy", "family", "uses_group_in_preprocessing", "description"),
  measure.vars = c(
    "median_run_partial_r2", "median_plate_partial_r2",
    "median_group_partial_r2", "median_absolute_sample_shift",
    "QC_median_feature_MAD", "QC_median_pairwise_spearman"
  ),
  variable.name = "metric", value.name = "value"
)
p_metrics <- ggplot(long_metrics, aes(strategy, value, fill = family)) +
  geom_col() +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  coord_flip() +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") +
  labs(title = "PACS strategy comparison", x = NULL, y = NULL)
ggsave(file.path(OUT_DIR, "plots", "strategy_metrics.png"), p_metrics, width = 12, height = 14, dpi = 180)

factor_plot_data <- factor_export
p_factor <- ggplot(factor_plot_data, aes(batch, k, fill = factor_type)) +
  geom_col(position = "dodge") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Estimated batch factors", x = NULL, y = "log2 offset")
ggsave(file.path(OUT_DIR, "plots", "batch_factors.png"), p_factor, width = 10, height = 6, dpi = 180)

best <- head(strategy_summary, 5)
report <- c(
  "# PACS multi-strategy benchmark",
  "",
  paste0("Generated: ", Sys.time()),
  "",
  "## Data used",
  "",
  paste0("- Study samples after README outlier removal: ", length(study_ids)),
  paste0("- Internal-QC samples: ", length(internal_ids)),
  paste0("- Plate-QC samples: ", length(plate_qc_ids)),
  paste0("- Protein groups: ", nrow(study_protein)),
  paste0("- Complete HELP set: ", length(complete_help)),
  paste0("- Missing-tolerant eligible HELP set: ", length(eligible_help)),
  "",
  "## Strategies",
  "",
  paste0("- ", strategy_info$strategy, ": ", strategy_info$description),
  "",
  "## Screening result",
  "",
  "The screening score is not a biological truth criterion. It rewards lower residual Run/Plate association, better QC behaviour, smaller unnecessary shifts, and preservation of the raw group effect.",
  "",
  paste0(
    seq_len(nrow(best)), ". **", best$strategy, "** — score ",
    signif(best$screening_score, 5),
    "; Run R2=", signif(best$median_run_partial_r2, 4),
    "; Plate R2=", signif(best$median_plate_partial_r2, 4),
    "; Group R2=", signif(best$median_group_partial_r2, 4)
  ),
  "",
  "## Interpretation rules",
  "",
  "- Do not select a method from PCA appearance alone.",
  "- Prefer methods that reduce Run/Plate effects in both study proteins and technical controls.",
  "- Treat group2-dependent preprocessing strategies as sensitivity analyses unless independently validated.",
  "- Review extreme sample factors and HELP detection counts before accepting sample-level correction.",
  "- The final choice should combine technical-control metrics, factor stability and biological plausibility."
)
writeLines(report, file.path(OUT_DIR, "REPORT.md"))

capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
message("PACS benchmark completed. Results written to: ", OUT_DIR)
