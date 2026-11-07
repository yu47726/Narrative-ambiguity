

# Narrative ambiguity and perspective-taking

Data and R code accompanying *Reduced narrative ambiguity amplifies the neural idiosyncrasy of mentalizing for sophisticated perspective-takers* by Yu Zhang, Bing Wu, and Yi Liu.

The shared files support statistical analyses of inter-subject dissimilarity in functional connectivity, strength centrality and eye-gaze trajectories, and inter-subject similarity of reactivation patterns. 

## Repository structure

All 14 input CSV files are directly inside `Data/`, and all four R scripts are directly inside `Scripts/`.

```text
Data/
  Study1_ISS_GlobalFC.csv
  Study1_ISS_GlobalSC.csv
  Study1_ISS_FCdiff.csv
  Study1_ISS_SCdiff.csv
  Study1_PT_Contr.csv
  Study1_EyeMove_ISD_mean.csv
  Study1_EyeMove_PT_Contr.csv
  Study2_ISS_GlobalFC.csv
  Study2_ISS_GlobalSC.csv
  Study2_ISS_FCdiff.csv
  Study2_ISS_SCdiff.csv
  Study2_PT_Contr.csv
  BESM_ISS_Network_Study1.csv
  BESM_ISS_Network_Study2.csv
Scripts/
  lme_utils.R
  run_global_LME.R
  run_connection_node_LME.R
  BESM_ISS_RunCompare.R
```

The scripts create `Results/` automatically; existing results are not required as input.

## Data

| File(s) in `Data/` | Contents and key columns |
|---|---|
| `Study1_ISS_GlobalFC.csv`, `Study2_ISS_GlobalFC.csv` | Global functional-connectivity dissimilarity: `sub1`, `sub2`, `Euc`, `Network`, `Ambiguity`. |
| `Study1_ISS_GlobalSC.csv`, `Study2_ISS_GlobalSC.csv` | Global strength-centrality dissimilarity, with the same columns. |
| `Study1_ISS_FCdiff.csv`, `Study2_ISS_FCdiff.csv` | Single-connection dissimilarity: `FCdiff` is the outcome; `FCName` identifies one of 15 connections. |
| `Study1_ISS_SCdiff.csv`, `Study2_ISS_SCdiff.csv` | Single-region strength-centrality dissimilarity: `Diff` is the outcome; `ROI` identifies one of six regions. |
| `Study1_EyeMove_ISD_mean.csv` | Eye-gaze trajectory dissimilarity: `sub1`, `sub2`, `ISD`, `Ambiguity`. |
| `Study1_PT_Contr.csv`, `Study1_EyeMove_PT_Contr.csv`, `Study2_PT_Contr.csv` | Participant scores: `pID`, perspective-taking (`IRI_PT`) and the IRI / Big Five subscales used as controls. |
| `BESM_ISS_Network_Study1.csv`, `BESM_ISS_Network_Study2.csv` | Reactivation-pattern similarity: `Dataset`, `sub1`, `sub2`, `all`, `upper`, `lower`, `Ambiguity`, `Net`. The model uses `lower` as the outcome and `upper` as a control. |

## Scripts

| File in `Scripts/` | Description |
|---|---|
| `lme_utils.R` | Shared functions for linear mixed effect model. |
| `run_global_LME.R` | Runs the global functional-connectivity, strength-centrality and eye-gaze dissimilarity analyses. |
| `run_connection_node_LME.R` | Fits separate ISD mixed-effects models for each connection or region; calculates and exports only the planned HH–LL interaction contrast. |
| `BESM_ISS_RunCompare.R` | Compares default-mode reactivation-pattern similarity between ambiguity conditions, controlling for upper-triangular BESM similarity. |
