## Global LME of the inter-subject dissimilarity (ISD): dyad group (HH / HL / LL) x ambiguity.

DataDir   <- "Data"
ResultDir <- file.path("Results", "ISS_LME", "Global")
CodeDir   <- "Scripts"
dir.create(ResultDir, recursive = TRUE, showWarnings = FALSE)

source(file.path(CodeDir, "lme_utils.R"))

datasets <- list(
  list(name = "Study1_GlobalFC", data_file = "Study1_ISS_GlobalFC.csv", scale_file = "Study1_PT_Contr.csv",
       dv_col = "Euc", network_val = "TOM"),
  list(name = "Study1_GlobalSC", data_file = "Study1_ISS_GlobalSC.csv", scale_file = "Study1_PT_Contr.csv",
       dv_col = "Euc", network_val = "TOM"),
  list(name = "Study1_EyeMove",  data_file = "Study1_EyeMove_ISD_mean.csv", scale_file = "Study1_EyeMove_PT_Contr.csv",
       dv_col = "ISD"),
  list(name = "Study2_GlobalFC", data_file = "Study2_ISS_GlobalFC.csv", scale_file = "Study2_PT_Contr.csv",
       dv_col = "Euc", network_val = "TOM"),
  list(name = "Study2_GlobalSC", data_file = "Study2_ISS_GlobalSC.csv", scale_file = "Study2_PT_Contr.csv",
       dv_col = "Euc", network_val = "TOM")
)

for (cfg in datasets) {
  cat("\n===", cfg$name, "===\n")
  res <- analyze_model(load_dyad_data(cfg, DataDir), cfg$dv_col)

  # FDR within each table of tests
  fdr_tables <- setdiff(names(res), c("ANOVA", "Params"))
  res[fdr_tables] <- lapply(res[fdr_tables], add_fdr)

  write.xlsx(res, file.path(ResultDir, paste0(cfg$name, "_global_LME.xlsx")))
}
cat("\nDone. Results saved to:", ResultDir, "\n")
