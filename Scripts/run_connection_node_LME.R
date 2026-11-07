## Connection / node-level LME of the ISD: the same model as run_global_LME.R, fitted separately to
## every functional connection (FCName) or region (ROI).
## One xlsx per dataset in Results/ISS_LME/Connection_Node. planned_contrast is the HH - LL planned contrast
## [HH - LL]Low-ambiguity - [HH - LL]High-ambiguity; FDR is applied over all connections (or regions).

DataDir   <- "Data"
ResultDir <- file.path("Results", "ISS_LME", "Connection_Node")
CodeDir   <- "Scripts"
dir.create(ResultDir, recursive = TRUE, showWarnings = FALSE)

source(file.path(CodeDir, "lme_utils.R"))

datasets <- list(
  list(name = "Study1_FCdiff", data_file = "Study1_ISS_FCdiff.csv", scale_file = "Study1_PT_Contr.csv",
       dv_col = "FCdiff", group_col = "FCName"),
  list(name = "Study1_SCdiff", data_file = "Study1_ISS_SCdiff.csv", scale_file = "Study1_PT_Contr.csv",
       dv_col = "Diff", group_col = "ROI"),
  list(name = "Study2_FCdiff", data_file = "Study2_ISS_FCdiff.csv", scale_file = "Study2_PT_Contr.csv",
       dv_col = "FCdiff", group_col = "FCName"),
  list(name = "Study2_SCdiff", data_file = "Study2_ISS_SCdiff.csv", scale_file = "Study2_PT_Contr.csv",
       dv_col = "Diff", group_col = "ROI")
)

for (cfg in datasets) {
  cat("\n===", cfg$name, "===\n")
  dyads <- load_dyad_data(cfg, DataDir)

  # one model per connection / region, tables stacked with a column identifying the connection / region
  per_group <- lapply(sort(unique(dyads[[cfg$group_col]])), function(grp) {
    cat("  -", grp, "\n")
    res <- analyze_model(dyads[dyads[[cfg$group_col]] == grp, ], cfg$dv_col, planned_only = TRUE)
    lapply(res, function(d) cbind(setNames(data.frame(grp), cfg$group_col), d))
  })
  res <- lapply(setNames(nm = names(per_group[[1]])), function(nm) do.call(rbind, lapply(per_group, `[[`, nm)))

  fdr_tables <- "planned_contrast"
  res[fdr_tables] <- lapply(res[fdr_tables], add_fdr)

  write.xlsx(res["planned_contrast"], file.path(ResultDir, paste0(cfg$name, "_LME.xlsx")))
}
cat("\nDone. Results saved to:", ResultDir, "\n")
