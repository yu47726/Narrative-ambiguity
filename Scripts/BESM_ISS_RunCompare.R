# Default-mode network only. Reporting direction: Low-ambiguity minus High-ambiguity; SE is not inflated.
rm(list = ls()) # clear

Idx <- "LowerCtrUp"

# Study 1: High- vs Low-ambiguity
file_name <- "Data/BESM_ISS_Network_Study1.csv"          # run from the repository root
SavePrefix <- "Results/Reactivation/LME_BESM_ISS_Network_Study1"

# Study 2: High- vs Low-ambiguity
# file_name <- "Data/BESM_ISS_Network_Study2.csv"
# SavePrefix <- "Results/Reactivation/LME_BESM_ISS_Network_Study2"

BESM_ISS_dyad <- read.csv(file_name)
BESM_ISS_dyad <- BESM_ISS_dyad[BESM_ISS_dyad$Net == "Default", , drop = FALSE]
stopifnot(nrow(BESM_ISS_dyad) > 0)
ambiguity_levels <- c("High-ambiguity", "Low-ambiguity")
stopifnot(setequal(unique(BESM_ISS_dyad$Ambiguity), ambiguity_levels))
BESM_ISS_dyad$Ambiguity <- factor(BESM_ISS_dyad$Ambiguity, levels = ambiguity_levels)
SaveFile <- paste0(SavePrefix, "_", Idx, ".RData")
SaveCsv <- paste0(SavePrefix, "_", Idx, ".csv")
dir.create(dirname(SavePrefix), recursive = TRUE, showWarnings = FALSE)


# Reuse dyad doubling only; keep this model and its original N-1 df convention.
source("Scripts/lme_utils.R")

stopifnot(all(c("sub1", "sub2", "Ambiguity", "Net", "upper", "lower") %in% colnames(BESM_ISS_dyad)))
n_unique_dyad <- nrow(unique(BESM_ISS_dyad[, c("sub1", "sub2")]))

BESM_ISS_dyad_double <- double_dyad_data(BESM_ISS_dyad)

# Set limits for df calculations
emm_options(pbkrtest.limit = 5000)
emm_options(lmerTest.limit = 5000)

dyad_BESM_ISS_RunCompare <- function(data)
{
  md <- lmer(scale(lower) ~ Ambiguity + scale(upper) + (1|sub1) + (1|sub2), data)

  
  emm_md <- emmeans(md, ~ Ambiguity)
  run_contrast <- contrast(emm_md, method = "revpairwise", adjust = "none")
  run_summary <- summary(run_contrast)
  
  B  <- run_summary$estimate[1]
  SE <- run_summary$SE[1]
  t  <- run_summary$t.ratio[1]
  p  <- run_summary$p.value[1]
  contrast_name <- as.character(run_summary$contrast[1])
  
  # Retain the original model-based SE; only the p-value df is adjusted.
  p_penalized  <- 2 * pt(-abs(t), (n_unique_dyad-1))
  
  emm_summary <- summary(emm_md)
  mean_Cond1 <- emm_summary$emmean[2]
  mean_Cond2 <- emm_summary$emmean[1]
  
  return(data.frame(
    Cond1 = ambiguity_levels[2],
    Cond2 = ambiguity_levels[1],
    Contrast = contrast_name,
    mean_Cond1 = mean_Cond1,
    mean_Cond2 = mean_Cond2,
    Cond1MinusCond2_B = B,
    Report_Contrast = "Low-ambiguity - High-ambiguity",
    Report_B = B,
    Report_SE = SE,
    Report_t = t,
    Report_df = n_unique_dyad - 1,
    Report_p = p_penalized,
    Cond1MinusCond2_t = t,
    Cond1MinusCond2_se = SE,
    Cond1MinusCond2_p = p,
    Cond1MinusCond2_p_penalized = p_penalized,
    mean_High = mean_Cond2,
    mean_Low = mean_Cond1,
    LowMinusHigh_B = B,
    LowMinusHigh_t = t,
    LowMinusHigh_se = SE,
    LowMinusHigh_p = p,
    LowMinusHigh_p_penalized = p_penalized
  ))
}

# Only Default remains after filtering; no network split is needed.
dyad_BESM_ISS_RunCompare_results <- data.frame(
  Net = "Default", dyad_BESM_ISS_RunCompare(BESM_ISS_dyad_double)
)
# Only Default is analyzed; BH across this single network equals the input p value.
dyad_BESM_ISS_RunCompare_results$Cond1MinusCond2_p_fdr <- p.adjust(
  dyad_BESM_ISS_RunCompare_results$Cond1MinusCond2_p,
  method = "BH"
)
dyad_BESM_ISS_RunCompare_results$Cond1MinusCond2_p_penalized_fdr <- p.adjust(
  dyad_BESM_ISS_RunCompare_results$Cond1MinusCond2_p_penalized,
  method = "BH"
)
# Preserve legacy aliases of the same contrast.
dyad_BESM_ISS_RunCompare_results$LowMinusHigh_p_fdr <-
  dyad_BESM_ISS_RunCompare_results$Cond1MinusCond2_p_fdr
dyad_BESM_ISS_RunCompare_results$LowMinusHigh_p_penalized_fdr <-
  dyad_BESM_ISS_RunCompare_results$Cond1MinusCond2_p_penalized_fdr
dyad_BESM_ISS_RunCompare_results$InputFile <- normalizePath(file_name, winslash = "/", mustWork = FALSE)

dyad_BESM_ISS_RunCompare_results

save(dyad_BESM_ISS_RunCompare_results, file = SaveFile)
write.csv(dyad_BESM_ISS_RunCompare_results, SaveCsv, row.names = FALSE)
