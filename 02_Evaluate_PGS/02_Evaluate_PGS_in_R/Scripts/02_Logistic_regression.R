# Run logistic regression

## ---- LoadDataFrame
load(file = "Data/Df.Rdata")
## ----

## ---- Exploration
# Data exploration before fitting regression models

# Check PGS before and after standardisation
Df %>%
  summarise(
    mean_PRS = mean(SCORESUM, na.rm = TRUE),
    sd_PRS = sd(SCORESUM, na.rm = TRUE),
    mean_PRS_standardised = mean(scale(SCORESUM), na.rm = TRUE),
    sd_PRS_standardised = sd(scale(SCORESUM), na.rm = TRUE)) %>% 
  kable()
  
  
# Plot distribution of standardised PGS split by epilepsy case/control
means <- Df %>%
  mutate(SCORESUM_scaled = scale(SCORESUM)) %>%
      group_by(Epilepsy) %>% 
      summarise(mean_value = mean(SCORESUM_scaled, na.rm = TRUE))
labels <- c("1" = "Case", "0" = "Control")
    
ggplot(Df, aes(x = scale(SCORESUM), colour = as.factor(Epilepsy), fill = as.factor(Epilepsy))) +
      geom_density(alpha = 0.6, position = "identity") +
      geom_vline(data = means, aes(xintercept = mean_value, colour = as.factor(Epilepsy)),
                 linetype = "dashed", linewidth = 1, alpha = 1,
                 show.legend = F) +
      scale_colour_viridis_d(name = "Epilepsy", labels = labels) +
      scale_fill_viridis_d(name = "Epilepsy", labels = labels) +
      scale_x_continuous(paste("Polygenic Score for epilepsy (standardised)")) +
      scale_y_continuous("") +
      theme_classic() +
      theme(text = element_text(family = "Calibri"),
            axis.text = element_text(size = 10),
            axis.title.x = element_text(size = 12, colour = "black", margin = margin(10,0,0,0)),
            axis.title.y = element_text(size = 12, colour = "black", margin = margin(0,10,0,0)),
            legend.title = element_blank(),
            legend.text = element_text(size = 10, colour = "black"),
            legend.position = "right")
## ----

## --- FitModels
logistic_regression_covar <- glm(Epilepsy ~ Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
                               family = binomial(link = 'logit'),
                               data = Df)

logistic_regression_PGS <- glm(Epilepsy ~ scale(SCORESUM) + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
                           family = binomial(link = 'logit'),
                           data = Df)

save(logistic_regression_covar, logistic_regression_PGS, file = "Results/logistic_regression_models.Rdata")
## ----

## ---- CheckModelCovar
performance::check_model(logistic_regression_covar)
simulateResiduals(logistic_regression_covar, plot=TRUE)
## ----

## ---- CheckModelFull
performance::check_model(logistic_regression_PGS)
simulateResiduals(logistic_regression_PGS, plot=TRUE)
## ----

## ---- RegressionResultsTable
tidy(logistic_regression_PGS, conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE) %>%
  filter(str_detect(term, "SCORESUM")) %>%
  kable(caption = "Results from logistic regression (back-transformation conducted so estimate = odds ratio")
## ----

## ---- RegressionResultsGraph
or_table <- tidy(logistic_regression_PGS, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95,) %>%
  filter(str_detect(term, "SCORESUM")) %>%
  mutate(
    term = str_replace(term, "scale\\(SCORESUM\\)", "PGS"),
    sig_label = ifelse(p.value < 0.05, "*", "")
  )

ggplot(or_table, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_y_continuous("Odds Ratio") +
  scale_x_discrete("Predictor") +
  # geom_text(aes(label = sig_label, y = conf.high + 0.05), size = 6) +  
  theme_classic() +
  theme(
    text = element_text(family = "Calibri"),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 12, colour = "black", margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(size = 12, colour = "black", margin = margin(0, 10, 0, 0)),
    legend.title = element_blank(),
    legend.text = element_text(size = 10, colour = "black"),
    legend.position = "right"
  )
  
## ----