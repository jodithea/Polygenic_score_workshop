# Calculate the area under the curve from receiver operating characteristic analysis
# Method 2: Run ROC analysis with covariates only and with covariates + PGS, then calculate AUC attributable to PGS

## ---- LoadDataFrame
load(file = "Data/Df.Rdata")
load(file = "Results/logistic_regression_models.Rdata")
## ----

## ---- ROCFitMethod2
# Use logistic regression models to run ROC analysis for covariates only and covariats + PGS
roc_model_covar <- roc(Df$Epilepsy, logistic_regression_covar$fitted.values)
roc_model_full <- roc(Df$Epilepsy, logistic_regression_PGS$fitted.values)

save(roc_model_covar, roc_model_full, file = "Results/roc_models_covar_full.Rdata")
## ----

## ---- AUCResultsMethod2
# DeLong's test for two correlated ROC curves
sig_test <- roc.test(roc_model_full, roc_model_covar)

tibble(
    AUC_covariates_only = as.numeric(auc(roc_model_covar)),
    AUC_covariates_plus_PGS = as.numeric(auc(roc_model_full)),
    AUC_change = AUC_covariates_only - AUC_covariates_plus_PGS,
    AUC_change_Z = sig_test[["statistic"]][["Z"]],
    AUC_change_p = sig_test[["p.value"]],
    AUC_change_CI_lower = sig_test$conf.int[1],
    AUC_change_CI_upper = sig_test$conf.int[2]) %>% 
    kable(caption = "AUC results and DeLong's test for two correlated ROC curves (does adding PGS to analysis significanlty increase AUC?)")
## ----

## ---- AUCPlotMethod2
roc_covar_df <- ggroc(roc_model_covar) %>% 
    .[["data"]] %>% 
    mutate(model = "covar")
  
roc_full_df  <- ggroc(roc_model_full) %>% 
    .[["data"]] %>% 
    mutate(model = "full")

roc_combined_df <- roc_covar_df %>% 
    bind_rows(roc_full_df)
  
auc_covar <- round(auc(roc_model_covar), 3)
auc_full  <- round(auc(roc_model_full), 3)
  
ggplot(roc_combined_df, aes(x = specificity, y = sensitivity, colour = model, linetype = model)) +
    geom_line() +
    scale_x_reverse("Specificity") +
    scale_y_continuous("Sensitivity") +
    scale_linetype_manual(values = c("covar" = "solid", "full" = "dashed")) +
    scale_colour_manual(values = c("covar" = "#160F3BFF", "full" = "#F4685CFF")) +
    annotate("text", 
             x = max(roc_combined_df$specificity, na.rm = TRUE) - 0.1, 
             y = min(roc_combined_df$sensitivity, na.rm = TRUE) + 0.05, hjust = 0, 
             label = paste("AUC (Covariates only) =", auc_covar), color = "#160F3BFF") +
    annotate("text", 
             x = max(roc_combined_df$specificity, na.rm = TRUE) - 0.1, 
             y = min(roc_combined_df$sensitivity, na.rm = TRUE) + 0.05 + 0.05, 
             hjust = 0, label = paste("AUC (Covariates + PGS) =", auc_full), color = "#F4685CFF") +
    theme_classic() +
    theme(text = element_text(family = "Calibri"),
          axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
          axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
          legend.title = element_blank(),
          legend.text = element_text(size = 10, colour = "black"),
          legend.position = "none")
  ## ----