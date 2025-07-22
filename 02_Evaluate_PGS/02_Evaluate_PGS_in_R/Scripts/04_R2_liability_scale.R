# Calculate variance explained in outcome variable risk on the liability scale attributable to PGS

## ---- LoadDataFrame
load(file = "Data/Df.Rdata")
## ----

## ---- VarExplainedLiabilityScaleFunction
# Variance explained by the PGS on the liability scale

# Convert observed R2 to the liability scale using poopulation prevalence of your disease/condition
# Use R code from Lee SH, et al., A Better Coefficient of Determination for Genetic Profile Analysis. Genetic Epidemiology, 2012. 36(3):214-224. 
# Made into a function
# lm = linear model
# K = population prevalence
# P = sample prevalence

mapToLiabilityScale = function(lm, K, P){
  thd = qnorm(1 - K) # the threshold on the normal distribution which truncates the proportion of disease prevalence K
  zv = dnorm(thd) # z (normal density)
  mv = zv/K # mean liability for case
  mv2 = -mv*K/(1-K) # mean liability for control
  
  y = lm$model[[1]]
  ncase = sum(y == 1)
  nt = length(y)
  ncont = nt - ncase
  R2O = var(lm$fitted.values)/(ncase/nt*ncont/nt) # R2 on the observed scale 
  
  theta = mv*(P-K)/(1-K)*(mv*(P-K)/(1-K)-thd) # theta in equation 15 of the publication
  cv = K*(1-K)/zv^2*K*(1-K)/(P*(1-P)) # C in equation 15 of the publication
  R2 = R2O*cv/(1+R2O*theta*cv)
  
  return(R2)
}
## ----

## ---- VarExplainedLiabilityScaleRun
# Calculate difference in R2 on liability scale of full model - covariates only model to identify variance attributable to PGS
# Bootstrap to get 90% CI
# Doing a one-sided test of significance (> 0) therefore use a one-sided 95% CI (equivalent to just the lower bound of a 90% CI)

set.seed(123)  # for reproducibility
n_boot <- 1000

r2_diffs <- numeric(n_boot)
    
for (i in 1:n_boot) {
      # Resample data with replacement
      boot_data <- Df[sample(1:nrow(Df), replace = TRUE), ]
      
      # Fit models with resampled data
      lm_covar <- lm(Epilepsy ~ Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
                   data = boot_data)

      lm_full <- lm(Epilepsy ~ SCORESUM + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
                    data = boot_data)
      
      # Calculate sample prevalence in the bootstrap sample
      ncase <- sum(boot_data$Epilepsy == 1)
      P_boot <- ncase / nrow(boot_data)
      
      # Liability R² (using lifetime population prevalence K as 0.012 (1.2%))
      r2_covar <- mapToLiabilityScale(lm_covar, K = 0.012, P = P_boot)
      r2_full <- mapToLiabilityScale(lm_full, K = 0.012, P = P_boot)
      
      r2_diffs[i] <- r2_full - r2_covar
    }
    
save(r2_diffs, file = "Results/r2_liability_attributable_PGS_bootstrapping.Rdata")
## ----

## ---- R2lResultsTable
# Doing a one-sided test of significance (> 0) therefore use a one-sided 95% CI (equivalent to just the lower bound of a 90% CI)
n_boot <- 1000
tibble(
  r2l_mean = mean(r2_diffs),
  r2l_SE = sd(r2_diffs) / sqrt(n_boot),
  r2l_Z = r2l_mean / r2l_SE,
  r2l_p = pnorm(r2l_Z, mean = 0, sd = 1, lower.tail = F),
  r2l_90perc_CI_lower = quantile(r2_diffs, 0.05),
  r2l_90perc_CI_upper = quantile(r2_diffs, 0.95)) %>% 
  kable(caption = "Variance explained in epilepsy risk on the liability scale attributable to epilepsy PGS and one-sided test of significance if R2 (liability) attributable to PGS is > 0")
## ----