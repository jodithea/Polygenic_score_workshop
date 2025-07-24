# Calculate the odds ratio of each decile compared to reference (lowest decile, or middle (5th and 6th) deciles)

## ---- LoadDataFrame
load(file = "Data/Df.Rdata")
## ----


## ---- TidyDataDecile
# Create deciles of PGS and set first decile as the reference/base group
Df_lowest <- Df %>%
  mutate(
    PGS_decile = ntile(scale(SCORESUM), 10),
    PGS_decile = factor(PGS_decile, levels = 1:10)
  )

# Create deciles of PGS and set middle deciles as the reference/base group
Df_mid <- Df %>%
  mutate(
    PGS_decile = ntile(scale(SCORESUM), 10),
    PGS_decile = case_when(
      PGS_decile %in% c("5", "6") ~ "5_6",
      TRUE ~ as.character(PGS_decile)
    ),
    PGS_decile = factor(PGS_decile, levels = c('5_6','1','2','3','4','7','8','9','10'))
  )
## ----

## ---- FitModelDecile
# Run logistic regression with the PGS decile as the predictor
decile_model_lowest <- glm(Epilepsy ~ PGS_decile + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
                    family = binomial(link = 'logit'),
                    data = Df_lowest)

decile_model_mid <- glm(Epilepsy ~ PGS_decile + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
                        family = binomial(link = 'logit'),
                        data = Df_mid)

save(decile_model_lowest, decile_model_mid, file = "Results/decile_model.Rdata")
## ----

## ---- CheckModelDecileLowest
performance::check_model(decile_model_lowest)
simulateResiduals(decile_model_lowest, plot=TRUE)
## ----

## ---- CheckModelDecileMid
performance::check_model(decile_model_mid)
simulateResiduals(decile_model_mid, plot=TRUE)
## ----

## ---- DecileModelResultsTableLowest
tidy(decile_model_lowest, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(str_detect(term, "PGS_decile")) %>%
  mutate(p.value_adjusted = p.adjust(p.value, method = "BH", n = 9)) %>% 
  kable(caption = "Results from logistic regression with decile of PGS as predictor and lowest decile as the reference (Back-transformation conducted so estimate = odds ratio. P-value adjusted for 9 comparisons with Benjamini-Hochberg method)")
## ----

## ---- DecileModelResultsTableMid
tidy(decile_model_mid, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(str_detect(term, "PGS_decile")) %>%
  mutate(p.value_adjusted = p.adjust(p.value, method = "BH", n = 8)) %>% 
  kable(caption = "Results from logistic regression with decile of PGS as predictor and middle deciles (5th and 6th) as the reference (Back-transformation conducted so estimate = odds ratio. P-value adjusted for 9 comparisons with Benjamini-Hochberg method)")
## ----

## ---- DecileModelResultsGraphLowest
or_table <- tidy(decile_model_lowest, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95,) %>%
  filter(str_detect(term, "PGS_decile")) %>%
  add_row(
    term = "PGS_decile1",
    estimate = 1,
    conf.low = 1,
    conf.high = 1,
    p.value = NA) %>% 
  mutate(
    decile = str_remove(term, "PGS_decile"),
    decile = factor(decile, levels = as.character(1:10)),
    p.value_adjusted = p.adjust(p.value, method = "BH", n = 9),
    sig_label = ifelse(p.value_adjusted < 0.05, "*", "")
  ) 

a <- ggplot(or_table, aes(x = decile, y = estimate)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_y_continuous("Odds Ratio") +
  scale_x_discrete("PGS Decile") +
  geom_text(aes(label = sig_label, y = conf.high + 1), size = 6) +  
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


## ---- DecileModelResultsGraphMid
or_table <- tidy(decile_model_mid, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95,) %>%
  filter(str_detect(term, "PGS_decile")) %>%
  add_row(
    term = "PGS_decile5",
    estimate = 1,
    conf.low = 1,
    conf.high = 1,
    p.value = NA) %>% 
  mutate(
    decile = str_remove(term, "PGS_decile"),
    decile = factor(decile, levels = as.character(1:10)),
    p.value_adjusted = p.adjust(p.value, method = "BH", n = 9),
    sig_label = ifelse(p.value_adjusted < 0.05, "*", "")
  ) 

a <- ggplot(or_table, aes(x = decile, y = estimate)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_y_continuous("Odds Ratio") +
  scale_x_discrete("PGS Decile",
                   labels = c('1', '2', '3', '4', '5 & 6', '7', '8', '9', '10')) +
  geom_text(aes(label = sig_label, y = conf.high + 1), size = 6) +  
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