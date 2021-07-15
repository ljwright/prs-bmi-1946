library(tidyverse)
library(haven)
library(glue)
library(broom)
library(gallimaufr)
library(magrittr)
library(ridittools)

rm(list = ls())

# 1. Load Data ----
# file_path <- "S:/IOECLS_Research/Bann/46c"
# 
# load_file <- function(file_name){
#   glue("{file_path}/{file_name}.dta") %>%
#     read_dta() %>%
#     rename_with(str_to_lower)
# }
# 
# set.seed(1)
# df_raw <- c("46c_ht_cog", 
#             "magels_fagels",
#             "bmi_khera2019_prs_withlabels_SCRAMBLED",
#             "2021_06_NSHD_two_childBMI_prs_scrambled_seqn340") %>%
#   map(load_file) %>%
#   reduce( ~ full_join(.x, .y, by = "nshdid_db1120", suffix = c("", "_xx")) %>%
#             select(-matches("_xx"))) %>%
#   mutate(id = sample(n()), .before = 1) %>%
#   select(-nshdid_db1120) %>%
#   arrange(id)
# save(df_raw, file = "Data/df_raw.Rdata")
load("Data/df_raw.Rdata")

# 2. Clean Data ----
year_to_age <- function(var_name){
  year <- str_sub(var_name, -2) %>% as.numeric()
  age <- ifelse(year < 46, 54 + year, year - 46)
  paste(str_sub(var_name, 1, -3), age, sep = "_")
}

sep_levels <- c("I Professional", "II Intermediate",
                "III Skilled Non-Manual", "III Skilled Manual",
                "IV", "V Unskilled")

df_long <- df_raw %>%
  rename(prs_k = bmi_prs_khera,
         prs_v = childbmi_prs_vogelezang2020, 
         prs_r = childbmi_prs_richardson2020) %>%
  rename_with(~ str_replace(.x, "n", ""), matches("^(h|w)tn")) %>%
  rename_with(~ str_replace(.x, "(u|x)$", ""), matches("^(ht|wt|bmi).(.*)$")) %>%
  rename_with(year_to_age, matches("^(ht|wt|bmi).(.*)$")) %>%
  mutate(across(matches("prs"), wtd_scale),
         across(matches("^(ht|wt|bmi)_"), as.numeric),
         across(matches("^(ht|wt|bmi)_"), 
                ~ case_when(.x < 0 ~ NA_real_,
                            .x >= 7777 ~ NA_real_,
                            TRUE ~ .x)),
         ht_63 = ht_63*100,
         survey_weight = as.numeric(inf),
         female = sex - 1,
         sep = as_factor(fsc50) %>% 
           fct_recode(NULL = "Unknown") %>%
           as.numeric() %>%
           sep_levels[.] %>%
           factor(levels = sep_levels),
         sep_ridit = toridit(table(sep))[sep] %>% as.numeric()) %>%
  select(id, survey_weight, matches("^(prs|bmi|wt|ht)_"), 
         female, sep, sep_ridit) %>%
  pivot_longer(matches("^(bmi|ht|wt)_"),
               names_to = c(".value", "age"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(age = as.numeric(age)) %>%
  arrange(id, age) %>%
  rename(height = ht, weight = wt) %>%
  zap_formats() %>%
  zap_label() %>%
  zap_labels()

# 3. BMI Height Correction ----
height_corr <- df_long %>%
  select(id, age, height, weight) %>%
  drop_na() %>%
  uncount(310, .id = "power") %>%
  mutate(power = (power-1)/100,
         bmi = weight / (height/100)^power) %>%
  group_by(age, power) %>%
  summarise(corr = cor(height, bmi),
            .groups = "drop")
  
height_power <- height_corr %>%
  group_by(age) %>%
  mutate(abs_corr = abs(corr)) %>%
  slice_min(abs_corr, n = 1) %>%
  ungroup() %>%
  select(age, power)

df_long <- df_long %>%
  left_join(height_power, by = "age") %>%
  mutate(bmi_corrected = weight / (height/100)^power) %>%
  select(-power)

save(df_long, file = "Data/df_long.Rdata")

# BMI Height Plots
height_corr %>%
  mutate(age = ordered(age)) %>%
  ggplot() +
  aes(x = power, y = corr) +
  facet_wrap(~ age) +
  geom_hline(yintercept = 0) +
  geom_line(data = rename(height_corr, age_f = age), 
            aes(group = age_f), color = "grey90") +
  geom_line(color = "#0072B2") +
  labs(x = "Power", y = "Correlation") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = FALSE)
ggsave("Images/height_power.png",
       width = 21, height = 21, units = "cm")
