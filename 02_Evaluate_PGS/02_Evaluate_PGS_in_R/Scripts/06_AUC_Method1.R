# Calculate the area under the curve from receiver operating characteristic analysis
# Method 1: Only include PGS in ROC model

## ---- LoadDataFrame
load(file = "Data/Df.Rdata")
## ----

## ---- ROCFitMethod1
model_PGS <- glm(Epilepsy ~ scale(SCORESUM), 
                 family = binomial(link = 'logit'), 
                 data = Df)

roc_model_PGS <- roc(Df$Epilepsy, model_PGS$fitted.values)
  
save(roc_model_PGS, file = "Results/roc_model_PGS.Rdata")
## ----

## ---- AUCResultsMethod1
# Doing a one-sided test of significance (> 0.5) therefore use a one-sided 95%CI (equivalent to just the lower bound of a 90% CI)
tibble(
      AUC = as.numeric(auc(roc_model_PGS)),
      AUC_SE = sqrt(var(roc_model_PGS)),
      AUC_Z = (AUC - 0.5) / AUC_SE,
      AUC_p = pnorm(AUC_Z, mean = 0, sd = 1, lower.tail = F),
      AUC_90perc_CI_lower = as.numeric(ci(roc_model_PGS, conf.level = 0.90))[1]) %>% 
  kable(caption = "AUC results and one-sided test of significance if AUC is > 0.5")
## ----

## ---- AUCPlotMethod1
plot(roc_model_PGS,
       print.auc = TRUE,
       main = paste("ROC curve for PGS of epilepsy"))
## ----