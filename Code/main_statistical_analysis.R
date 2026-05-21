# ==============================================================================
# Title: Comparative Effectiveness of MET-TKIs vs Systemic Therapies in METex14 NSCLC
# Description: Core statistical pipeline for pseudo-IPD reconstruction, covariate 
#              balancing (IPTW/Overlap Weighting), RMST, and sensitivity analysis.
# Authors: Fangfang Shen, Zheqing Zhu, et al.
# Dependencies: dplyr, readr, survival, cobalt, EValue
# Input Files: final_combined_analysis_data_fixed.csv, master_baseline_table_fixed.csv
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(cobalt)
  library(EValue)
})

# Set working directory and load inputs
# ⚠️ USER NOTE: Please update the path below to your local data directory before running.
setwd("/Data") 
cat("\n[INIT] Loading pseudo-IPD and aggregate baseline matrices...\n")

df_ipd <- read_csv("final_combined_analysis_data_fixed.csv", show_col_types = FALSE) 
df_base <- read_csv("master_baseline_table_fixed.csv", show_col_types = FALSE)                                  

# ==============================================================================
# 1. Helper Function: Pseudo-IPD Simulation & Propensity Score Weighting
# ==============================================================================
# This function performs a dual-key merge to prevent Cartesian duplication across 
# survival endpoints. It applies propensity score weighting based on user-defined 
# estimands (ATT via IPTW, or ATO via Overlap Weighting).
generate_weighted_cohort <- function(df_ipd, df_base, target_endpoint, method = "IPTW", exclude_024 = FALSE) {
  
  # Isolate specific endpoint and define control pool
  df_sub <- df_ipd %>% filter(endpoint == target_endpoint)
  if (exclude_024) {
    df_sub <- df_sub %>% filter(cohort == "Experimental" | (cohort == "Control" & !grepl("KEYNOTE-024", group)))
  } else {
    df_sub <- df_sub %>% filter(cohort == "Experimental" | grepl("KEYNOTE-024", group))
  }
  
  # Merge baseline characteristics using dual keys (Treatment + Endpoint)
  # Missing baseline proportions are imputed with neutral/median values to maintain cohort sizes
  df_full <- df_sub %>% 
    left_join(df_base, by = c("group" = "Treatment", "endpoint" = "Endpoint")) %>%
    mutate(
      Median_Age = ifelse(is.na(Median_Age), median(Median_Age, na.rm = TRUE), Median_Age),
      Male_pct = ifelse(is.na(Male_pct), 50, Male_pct),
      # 【核心对齐 1】：匹配定稿论文的插补比例，确保生成一致的虚拟队列
      ECOG_ge_1_pct = ifelse(is.na(ECOG_ge_1_pct), 70, ECOG_ge_1_pct),
      Never_Smoker_pct = ifelse(is.na(Never_Smoker_pct), 40, Never_Smoker_pct),
      Brain_Mets_pct = ifelse(is.na(Brain_Mets_pct), 15, Brain_Mets_pct)
    )
  
  # Reconstruct baseline covariates for Propensity Score estimation
  # Seed is fixed to guarantee exact reproducibility of the Monte Carlo simulation
  set.seed(2026) 
  df_sim <- df_full %>% rowwise() %>%
    mutate(
      # 【核心对齐 2】：严格匹配定稿论文的 RNG 消耗顺序 (先 rbinom, 后 rnorm)
      Sim_Male = rbinom(1, 1, Male_pct/100),
      Sim_ECOG = rbinom(1, 1, ECOG_ge_1_pct/100),
      Sim_Smoker = rbinom(1, 1, Never_Smoker_pct/100),
      Sim_BrainMets = rbinom(1, 1, Brain_Mets_pct/100),
      Sim_Age = round(rnorm(1, mean = Median_Age, sd = 8), 1)
    ) %>%
    ungroup() %>%
    mutate(Treatment_Binary = ifelse(cohort == "Experimental", 1, 0))
  
  # Logistic regression for PS estimation (na.action added for strict robustness)
  ps_mod <- glm(Treatment_Binary ~ Sim_Age + Sim_Male + Sim_ECOG + Sim_Smoker + Sim_BrainMets, 
                family = binomial(link = "logit"), data = df_sim, na.action = na.exclude)
  
  df_sim$PS <- pmax(0.01, pmin(0.99, predict(ps_mod, type = "response")))
  
  # Compute respective weights based on causal estimand
  if (method == "IPTW") {
    # Average Treatment Effect on the Treated (ATT)
    # Includes 99th percentile truncation to mitigate variance inflation from extreme weights
    df_sim$Weight <- ifelse(df_sim$Treatment_Binary == 1, 1, df_sim$PS / (1 - df_sim$PS))
    df_sim$Weight <- pmin(df_sim$Weight, quantile(df_sim$Weight, 0.99))
  } else if (method == "OW") {
    # Average Treatment Effect on the Overlap population (ATO)
    df_sim$Weight <- ifelse(df_sim$Treatment_Binary == 1, 1 - df_sim$PS, df_sim$PS)
  }
  
  return(df_sim)
}

# ==============================================================================
# 2. Frontline Efficacy: MET-TKIs vs. Standard Platinum-Based Chemotherapy
# ==============================================================================
cat("\n--- Section 2: 1L PFS Efficacy vs. Standard Chemotherapy (IPTW/ATT) ---\n")
df_1L_pfs <- generate_weighted_cohort(df_ipd, df_base, "PFS1", method = "IPTW", exclude_024 = TRUE)

# Diagnostic: Verify covariate balance (Target SMD < 0.20 for all variables)
bal_obj_1L <- bal.tab(Treatment_Binary ~ Sim_Age + Sim_Male + Sim_ECOG + Sim_Smoker + Sim_BrainMets, 
                      data = df_1L_pfs, weights = df_1L_pfs$Weight, method = "weighting", 
                      estimand = "ATT", un = TRUE, continuous = "std")
print(bal_obj_1L$Balance[, c("Diff.Un", "Diff.Adj")])

# Fit weighted Cox proportional hazards model using robust standard errors
fit_cox_1L <- coxph(Surv(time, status) ~ Treatment_Binary, data = df_1L_pfs, weights = Weight, robust = TRUE)
sum_cox_1L <- summary(fit_cox_1L)
cat(sprintf("Adjusted HR for 1L PFS: %.2f (95%% CI: %.2f - %.2f), Robust P = %.3e\n", 
            sum_cox_1L$conf.int[1, "exp(coef)"], sum_cox_1L$conf.int[1, "lower .95"], 
            sum_cox_1L$conf.int[1, "upper .95"], sum_cox_1L$coefficients[1, "Pr(>|z|)"]))

# ==============================================================================
# 3. Immunotherapy Comparison: Overlap Weighting and RMST Dynamics
# ==============================================================================
cat("\n--- Section 3: 1L PFS Efficacy vs. Pembrolizumab (OW/ATO) ---\n")
df_024_pfs <- generate_weighted_cohort(df_ipd, df_base, "PFS1", method = "OW", exclude_024 = FALSE)

# Impose an administrative censoring threshold at 36 months to align follow-up windows
cutoff_months <- 36
df_024_pfs <- df_024_pfs %>%
  mutate(
    status = ifelse(time > cutoff_months, 0, status),
    time = ifelse(time > cutoff_months, cutoff_months, time)
  )

# Diagnostic: Assess exact covariate balance under OW framework (Target SMD ~ 0.000)
bal_obj_OW <- bal.tab(Treatment_Binary ~ Sim_Age + Sim_Male + Sim_ECOG + Sim_Smoker + Sim_BrainMets, 
                      data = df_024_pfs, weights = df_024_pfs$Weight, method = "weighting", 
                      estimand = "ATO", un = TRUE, continuous = "std")
print(bal_obj_OW$Balance[, c("Diff.Un", "Diff.Adj")])

# Calculate Restricted Mean Survival Time (RMST) to address non-proportional hazards
calculate_weighted_rmst <- function(df, tau) {
  fit <- survfit(Surv(time, status) ~ Treatment_Binary, data = df, weights = Weight)
  tbl <- summary(fit, rmean = tau)$table
  
  rmst_0 <- tbl["Treatment_Binary=0", "rmean"]
  se_0   <- tbl["Treatment_Binary=0", "se(rmean)"]
  rmst_1 <- tbl["Treatment_Binary=1", "rmean"]
  se_1   <- tbl["Treatment_Binary=1", "se(rmean)"]
  
  # Wald test for RMST difference incorporating robust SEs
  diff <- rmst_1 - rmst_0
  se_diff <- sqrt(se_0^2 + se_1^2)
  ci_l <- diff - 1.96 * se_diff
  ci_u <- diff + 1.96 * se_diff
  pval <- 2 * pnorm(-abs(diff / se_diff))
  
  return(c(Diff = diff, Lower = ci_l, Upper = ci_u, Pval = pval))
}

cat("Evaluating temporal differences in PFS using RMST:\n")
for (m in c(12, 24, 36)) {
  res <- calculate_weighted_rmst(df_024_pfs, m)
  cat(sprintf("Milestone %d Months | Difference: %+.2f (95%% CI: %.2f to %.2f), P = %.3f\n", 
              m, res["Diff"], res["Lower"], res["Upper"], res["Pval"]))
}

# ==============================================================================
# 4. Previously Treated Efficacy: Salvage Setting (2L+)
# ==============================================================================
cat("\n--- Section 4: 2L+ Overall Survival vs. Salvage Regimens ---\n")
df_2L_os <- generate_weighted_cohort(df_ipd, df_base, "OS2", method = "IPTW", exclude_024 = TRUE)

fit_cox_2L <- coxph(Surv(time, status) ~ Treatment_Binary, data = df_2L_os, weights = Weight, robust = TRUE)
sum_cox_2L <- summary(fit_cox_2L)
cat(sprintf("Adjusted HR for 2L+ OS: %.2f (95%% CI: %.2f - %.2f), Robust P = %.3e\n", 
            sum_cox_2L$conf.int[1, "exp(coef)"], sum_cox_2L$conf.int[1, "lower .95"], 
            sum_cox_2L$conf.int[1, "upper .95"], sum_cox_2L$coefficients[1, "Pr(>|z|)"]))

# ==============================================================================
# 5. Quantitative Bias Analysis: E-value for Unmeasured Confounding
# ==============================================================================
cat("\n--- Section 5: Robustness Assessment (E-value) ---\n")

gold_standard_hr <- sum_cox_1L$conf.int[1, "exp(coef)"]
gold_standard_lower <- sum_cox_1L$conf.int[1, "lower .95"]
gold_standard_upper <- sum_cox_1L$conf.int[1, "upper .95"]

e_val_res <- evalues.HR(est = gold_standard_hr, lo = gold_standard_lower, hi = gold_standard_upper, rare = FALSE)

# Statistical Note: For protective treatments (HR < 1), the closest limit to the null (HR=1) 
# is the upper confidence limit. Thus, the E-value's lower bound corresponds to the HR's upper limit.
e_value_estimate <- e_val_res[2, "point"]
e_value_lower_bound <- e_val_res[2, "upper"] 

cat(sprintf("E-value Point Estimate: %.2f\nLower Confidence Limit: %.2f\n", 
            e_value_estimate, e_value_lower_bound))
cat("\n[COMPLETE] Pipeline executed successfully.\n")
