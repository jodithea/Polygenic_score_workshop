## ---- LoadPackages
library(knitr)        # for kable tables
library(performance)  # for model checking
library(DHARMa)       # for model checking
library(effects)      # for nice summary of model results (allEffects)
library(broom.mixed)  # for tidy output
library(ggplotify)    # to use patchwork on lattice plots (created by allEffects)
library(MuMIn)        # to calculate R2
library(pROC)         # for ROC curves and AUC
library(patchwork)    # for combining ggplot plots
library(tidyverse)
## ----

## ---- LoadData
# Update this when have final results from workshop:
PGS_df <- read.table(file = "/path/PGS_Workshop/01_Create_PGS/Plink_PGS_scores/ILAE3_Caucasian_all_epilepsy_EUR.QC_PGS_scores.profile", header = T)

# Read in covariates data
cov_df <- read.table(file = "/path/PGS_Workshop/Genotype_data/EUR.cov", header = T)

# Read in ancestry PCs data
PCs_df <- read.table(file = "/path/PGS_Workshop/Genotype_data/EUR.eigenvec", header = F, col.names = c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6"))
  
# Relatedness data
related_df <- read.table(file = "/path/PGS_Workshop/02_Evaluate_PGS/EUR_relatedness.genome", header = T)
## ----

## ---- TidyData
# Merge these three dataframes (PGS, covariates and PCs) and add simulated data for our phenotype = epilepsy (case/control)

set.seed(123) # set seed for reproducibility
Df <- PGS_df %>% 
  left_join(cov_df) %>% 
  left_join(PCs_df) %>% 
  mutate(
    lp = -1.5 + 0.8 * scale(SCORESUM) + 
      0.05 * Sex + 0.01 * PC1 + 0.005 * PC2 +
      0.005 * PC3 + 0.00001 * PC4 + 0.0003 * PC5 + 0.00008 * PC6,
    prob = plogis(lp),
    Epilepsy = rbinom(n(), size = 1, prob = prob)) %>% 
  select(-c(lp, prob))



# Remove 1 individual from each pair that have pi_hat > 0.1875 with controls preferentially removed.
# This removes individuals that are first- and second-degress relatives (0.1875 is halfway between 2nd and 3rd degree relatives)
# Note in our data that related_df is empty - no individuals have pi_hat > 0.09 (this dataset was QCed and related individuals were already removed)
# So this code is commented out, but you can use it if your data includes related individuals

# # Df with case/ctrl status for each IID
# Epilepsy_status <- Df %>% 
#   select(IID, Epilepsy)
# 
# # Join related df with CP Epilepsy_status so have case/control status for both IID1 and IID2
# 
# related_epilepsy_status <- related_df %>%
#   left_join(Epilepsy_status, by = c("IID1" = "IID")) %>%
#   rename(Epilepsy_IID1 = Epilepsy) %>%
#   left_join(Epilepsy_status, by = c("IID2" = "IID")) %>%
#   rename(Epilepsy_IID2 = Epilepsy)
# 
# # Filter pairs with PI_HAT > 0.1875, preferentially removing controls and keeping cases
# related_remove <- related_epilepsy_status %>%
#   filter(PI_HAT > 0.1875) %>%
#   rowwise() %>% 
#   mutate(
#     remove = case_when(
#       Epilepsy_IID1 == 0 & Epilepsy_IID2 == 1 ~ IID1,  # Prefer removing control (IID1)
#       Epilepsy_IID1 == 1 & Epilepsy_IID2 == 0 ~ IID2,  # Prefer removing control (IID2)
#       Epilepsy_IID1 == 0 & Epilepsy_IID2 == 0 ~ sample(c(IID1, IID2), 1),  # Randomly remove one if both are controls
#       Epilepsy_IID1 == 1 & Epilepsy_IID2 == 1 ~ sample(c(IID1, IID2), 1)   # Randomly remove one if both are cases
#     )
#   ) %>% 
#   select(remove)
# 
# # Remove these individuals from df
# Df <- Df %>% 
#   anti_join(related_remove, by = join_by("IID" == "remove"))

save(Df, file = "Data/Df.Rdata")
## ----

