#!/usr/bin/env Rscript

library(tidyverse)
library(broom)
library(glue)

# --- Configuration ---
PREDICTIONS_FILE <- "../Data/Organ_Age_Predictions.csv"
METADATA_FILE    <- "../Data/Metadata.csv"
OUTPUT_FILE      <- "../Results/agegap_sex_adjusted.csv"

TARGET_ORGANS <- c("Immune", "Brain", "Liver", "Artery", "Adipose",
                   "Pancreas", "Heart", "Kidney", "Lung", "Muscle", "Intestine")


# --- Main Script ---

# Load and prepare input data
predictions <- read_csv(PREDICTIONS_FILE, show_col_types = FALSE)

# Load metadata with Sex
sex_table <- read_csv(METADATA_FILE, show_col_types = FALSE) %>%
  distinct(ParticipantID, Sex)

# Join predictions with sex information and filter for target organs
predictions_with_sex <- predictions %>%
  filter(Organ %in% TARGET_ORGANS) %>%
  left_join(sex_table, by = "ParticipantID") %>%
  drop_na(Sex)

# Compute LOESS residuals to define age gaps
# Residuals are calculated per organ, cross-validation fold, and visit round.
loess_residuals <- predictions_with_sex %>%
  group_by(Organ, Fold, Round) %>%
  filter(n() > 2) %>%
  reframe({
    fit <- loess(PredictedAge ~ Age,
                 data    = cur_data_all(),
                 span    = 1.0,
                 control = loess.control(surface = "direct"))
    broom::augment(fit, newdata = cur_data_all())
  }) %>%
  rename(AgeGap = .resid) %>%
  ungroup()

# Z-score the age gaps within each organ and visit round
zscored_age_gaps <- loess_residuals %>%
  group_by(Organ, Round) %>%
  mutate(AgeGap_Z = as.numeric(scale(AgeGap))) %>%
  ungroup()

# Reshape data to a wide format for baseline (R1) and follow-up (R2)
zscored_gaps_wide <- zscored_age_gaps %>%
  filter(Round %in% 1:2) %>%
  group_by(ParticipantID, Organ) %>%
  filter(n_distinct(Round) == 2) %>% # Ensure pairs of visits are present
  summarise(
    g_R1   = AgeGap_Z[Round == 1],
    g_R2   = AgeGap_Z[Round == 2],
    Sex    = first(Sex),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols     = c(ParticipantID, Sex),
    names_from  = Organ,
    values_from = c(g_R1, g_R2),
    names_glue  = "{Organ}_{.value}"
  ) %>%
  drop_na()

# Adjust baseline (R1) z-scored age gaps for sex
baseline_cols <- grep("_g_R1$", names(zscored_gaps_wide), value = TRUE)
baseline_residuals <- map_dfc(baseline_cols, function(col) {
  X     <- model.matrix(~ Sex, data = zscored_gaps_wide)
  resid <- zscored_gaps_wide[[col]] - X %*% coef(lm.fit(X, zscored_gaps_wide[[col]]))
  tibble(!!glue("{col}_res") := as.numeric(resid))
})

# Adjust follow-up (R2) z-scored age gaps for sex
followup_cols <- grep("_g_R2$", names(zscored_gaps_wide), value = TRUE)
followup_residuals <- map_dfc(followup_cols, function(col) {
  X     <- model.matrix(~ Sex, data = zscored_gaps_wide)
  resid <- zscored_gaps_wide[[col]] - X %*% coef(lm.fit(X, zscored_gaps_wide[[col]]))
  tibble(!!glue("{col}_res") := as.numeric(resid))
})

# Combine identifiers with the new sex-adjusted residuals
final_adjusted_gaps <- bind_cols(
  zscored_gaps_wide["ParticipantID"],
  baseline_residuals,
  followup_residuals
)

# Save the final results to a file
write_csv(final_adjusted_gaps, OUTPUT_FILE)

print(glue("Sex-adjusted age gaps saved to {OUTPUT_FILE}"))