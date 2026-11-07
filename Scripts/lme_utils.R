## lme_utils.R
## Shared functions of the ISD LME analyses (sourced by run_global_LME.R and run_connection_node_LME.R).
##
## Model: scale(ISD) ~ Ambiguity * (PT_dyad + covariates) + (1|sub1) + (1|sub2), fitted on the doubled data.
##
## DF-only adjustment for the doubled data (SE is retained from the fitted model):
##   N = number of unique dyads, k = rank of the fixed-effects design matrix (intercept included)
##   df  : N-k for the denominator of every F test and for every contrast
##   SE  : original model-based SE; no multiplicative adjustment
##   p   : computed from the raw t / F with N-k df
##   FDR : Benjamini-Hochberg on the corrected p (add_fdr)
## One N-k serves all tests: Ambiguity (within dyads) has only 2 levels, so the error df of the within
## effects (Ambiguity, Ambiguity x PT), (a-1)(N-g) with a = 2, equals that of the PT effect, N-g.

library(lmerTest)
library(emmeans)
library(openxlsx)

## Preserve values and row order; append the dyads with participant roles swapped.
double_dyad_data <- function(data) {
  stopifnot(all(c("sub1", "sub2") %in% names(data)))
  rbind(data, transform(data, sub1 = sub2, sub2 = sub1))
}

## Reads one dataset (cfg: data_file, scale_file, optional network_val) and returns the doubled dyad data.
## For every scale variable in the participant-level file (all columns except pID), the participants are
## split at the median (strictly greater than the median = "high") and every dyad is assigned to
## "high_high" (HH), "low_low" (LL) or "low_high" (mixed, HL): factor column dyad_<scale>_binary.

load_dyad_data <- function(cfg, data_dir) {
  data_df  <- read.csv(file.path(data_dir, cfg$data_file))
  scale_df <- read.csv(file.path(data_dir, cfg$scale_file))
  if (!is.null(cfg$network_val)) data_df <- data_df[data_df$Network == cfg$network_val, ]
  data_df$Ambiguity <- as.factor(data_df$Ambiguity)
  data_df <- data_df[data_df$sub1 %in% scale_df$pID & data_df$sub2 %in% scale_df$pID, ]   # both members need scores

  for (var in setdiff(names(scale_df), "pID")) {
    high   <- setNames(scale_df[[var]] > median(scale_df[[var]], na.rm = TRUE), scale_df$pID)
    n_high <- high[data_df$sub1] + high[data_df$sub2]
    data_df[[paste0("dyad_", var, "_binary")]] <- factor(c("low_low", "low_high", "high_high")[n_high + 1])
  }

  double_dyad_data(data_df)
}

## Correction parameters of one fitted model: list(N, k, df)
chen_params <- function(model, data) {
  a <- as.character(data$sub1)
  b <- as.character(data$sub2)
  N <- length(unique(ifelse(a < b, paste(a, b), paste(b, a))))   # unordered sub1/sub2 pairs
  k <- qr(getME(model, "X"))$rank
  list(N = N, k = k, df = N - k)
}

## emmeans contrast table
correct_contrast <- function(d, chen) {
  d$t.ratio <- d$estimate / d$SE                  # t and CI both use the original model-based SE
  d$z.ratio <- NULL                               # emmeans reports z.ratio because the df is Inf here
  d$df      <- chen$df
  d$p.value <- 2 * pt(-abs(d$t.ratio), chen$df)
  d$CI_low  <- d$estimate - qt(0.975, chen$df) * d$SE
  d$CI_high <- d$estimate + qt(0.975, chen$df) * d$SE
  d
}

add_fdr <- function(d) {
  d$p.fdr <- p.adjust(d$p.value, method = "fdr")
  d
}

## Fits the LME to the doubled data and returns all corrected result tables:
##   ANOVA, SimpleEffect_PT_by_Ambiguity (F test of PT within each level of Ambiguity),
##   PT_in_Ambiguity (HH-HL, HH-LL, HL-LL within each level), Ambiguity_in_PT (Low - High within each group),
##   planned_contrast ([HH - LL]Low-ambiguity - [HH - LL]High-ambiguity)
analyze_model <- function(data, dv_col, iv_col = "dyad_IRI_PT_binary", planned_only = FALSE) {
  ctr_cols <- setdiff(grep("^dyad_.*_binary$", names(data), value = TRUE), iv_col)
  formula  <- paste0("scale(", dv_col, ") ~ Ambiguity * (", paste(c(iv_col, ctr_cols), collapse = " + "),
                     ") + (1|sub1) + (1|sub2)")
  model    <- lmer(as.formula(formula), data = data, REML = FALSE)
  chen     <- chen_params(model, data)

  if (planned_only) {
    # Retain the full fitted model; evaluate only the HH and LL marginal means.
    emm <- emmeans(model, as.formula(paste("~ Ambiguity *", iv_col)),
                   at = setNames(list(c("high_high", "low_low")), iv_col),
                   lmer.df = "asymptotic")
    planned <- correct_contrast(as.data.frame(contrast(
      emm, interaction = c("revpairwise", "pairwise"), adjust = "none"
    )), chen)
    return(list(planned_contrast = planned))
  }

  anova_df <- as.data.frame(anova(model))                       # Type III, Satterthwaite
  anova_df <- data.frame(Effect = rownames(anova_df), anova_df, check.names = FALSE, row.names = NULL)
  anova_df$DenDF       <- chen$df
  anova_df[["Pr(>F)"]] <- pf(anova_df[["F value"]], anova_df$NumDF, chen$df, lower.tail = FALSE)

  # df is replaced by N-k below, so the cheap asymptotic df is enough here
  emm <- emmeans(model, as.formula(paste("~ Ambiguity *", iv_col)), lmer.df = "asymptotic")
  contrasts <- function(...) correct_contrast(as.data.frame(contrast(emm, ..., adjust = "none")), chen)

  simple <- as.data.frame(joint_tests(emm, by = "Ambiguity"))
  names(simple)[names(simple) == "model term"] <- "Effect"
  simple$df2     <- chen$df
  simple$p.value <- pf(simple$F.ratio, simple$df1, chen$df, lower.tail = FALSE)

  # Keep only the prespecified HH - LL difference-in-differences.
  planned <- contrasts(interaction = c("revpairwise", "pairwise"))
  planned <- planned[planned[[paste0(iv_col, "_pairwise")]] == "high_high - low_low", , drop = FALSE]

  list(ANOVA                        = anova_df,
       SimpleEffect_PT_by_Ambiguity = simple[c("Ambiguity", "Effect", "df1", "df2", "F.ratio", "p.value")],
       PT_in_Ambiguity              = contrasts("pairwise", by = "Ambiguity"),
       Ambiguity_in_PT              = contrasts("revpairwise", by = iv_col),
       planned_contrast             = planned,
       Params                  = as.data.frame(chen))
}
