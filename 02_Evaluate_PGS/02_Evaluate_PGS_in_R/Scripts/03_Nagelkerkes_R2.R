# Calculate Nagelkerke's R2 for variance explained by PGS

## ---- LoadDataFrameAndModels
load(file = "Data/Df.Rdata")
load(file = "Results/logistic_regression_models.Rdata")
## ----

## ---- NagelkerkeR2
r2_PGS <- r2_nagelkerke(logistic_regression_PGS) - r2_nagelkerke(logistic_regression_covar)

r2_PGS %>% 
  kable(caption = "Nagelkerke's R2 attributable to the PGS")
## ----