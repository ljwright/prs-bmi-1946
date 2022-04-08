library(tidyverse)
library(haven)
library(glue)
library(broom)
library(gallimaufr)
library(magrittr)
library(ridittools)
library(scales)
library(sjmisc)
library(corrr)

rm(list = ls())

# 1. Load Data ----
file_path <- Sys.getenv("raw_fld")

load_file <- function(file_name){
  glue("{file_path}/{file_name}.dta") %>%
    read_dta() %>%
    rename_with(str_to_lower)
}

set.seed(1)
file_names <- c("survey_data", 
                "bmi_khera2019_prs_withlabels_SCRAMBLED",
                "2022_02_NSHD_adultBMI_prs_Richardson2020_scrambled_seqn340_id",
                "2021_06_NSHD_two_childBMI_prs_scrambled_seqn340",
                "2022_01_NSHD_BMI_data_for_D_Bann_inc_Plink_format_IID_and_10_PCs_SCRAMBLED")

df_raw <- map(file_names, load_file) %>%
  reduce( ~ full_join(.x, .y, by = "nshdid_db1120", suffix = c("", "_xx")) %>%
            select(-matches("_xx"))) %>%
  mutate(id = sample(n()), .before = 1) %>%
  select(-nshdid_db1120) %>%
  arrange(id)
save(df_raw, file = "Data/df_raw.Rdata")

rm(file_names, file_path)

# 2. Clean Data ----
load("Data/df_raw.Rdata")

# Variable Labels
do_file <- read_lines("Data/STATA_Script_dbannZZdhxddp_20220216-115754.do") %>%
  tibble(line = .)

value_labels <- do_file %>%
  filter(str_detect(line, "^label define")) %>%
  mutate(line = str_replace(line, "label define ", "") %>%
           str_replace(", modify", "") %>%
           str_replace_all('\\"', "XX")) %>%
  separate(line, c("line", "label"), sep = " XX", extra = "drop") %>%
  mutate(label = str_replace(label, "XX", "")) %>%
  separate(line, c("var", "value"), sep = " ") %>%
  mutate(value = str_replace_all(value, "\\'", "") %>%
           as.integer())

variable_labels <- do_file %>%
  filter(str_detect(line, "^label variable")) %>%
  mutate(line = str_replace(line, "label variable ", "") %>%
           str_replace_all('\\"', "XX")) %>%
  separate(line, c("variable", "label"), sep = " XX", extra = "drop") %>%
  mutate(label = str_replace(label, "XX", ""))

# Functions
year_to_age <- function(var_name){
  year <- str_sub(var_name, -2) %>% as.numeric()
  age <- ifelse(year < 46, 54 + year, year - 46)
  paste(str_sub(var_name, 1, -3), age, sep = "_")
}

sep_levels <- c("I Professional", "II Intermediate",
                "III Skilled Non-Manual", "III Skilled Manual",
                "IV Semi-Skilled", "V Unskilled")

educ_levels <- c("No Qualifications" = 0, "Vocational Course, Proficiency Only" = 1,
                 "Below GCE or Burnham C" = 2, "O Level or Burnham B" = 3,
                 "A Level or Burnham B" = 4, "Burnham A2" = 5, "Degree or Equivalent" = 6,
                 "Masters Degree" = 7, "Doctorate" = 8)

make_ridit <- function(x, rev = FALSE){
  fac <- factor(x)
  if (rev) fac <- fct_rev(fac)
  toridit(table(fac))[fac] %>%
    as.numeric()
}

clean_between <- function(x, low, high) ifelse(between(x, low, high), x, NA)
truncate <- function(x, a, b) pmin(pmax(x, a), b)

# Clean
df_long <- df_raw %>%
  # Polygenic Scores
  rename(prs_r_adult  = adultbmi_prs_richardson2020,
         prs_k = bmi_prs_khera,
         prs_v = childbmi_prs_vogelezang2020, 
         prs_r_child = childbmi_prs_richardson2020) %>%
  mutate(across(matches("prs"), wtd_scale)) %>%
  # Anthropometrics
  rename_with(~ str_replace(.x, "n", ""), matches("^(h|w)tn")) %>%
  rename_with(~ str_replace(.x, "(u|x)$", ""), matches("^(ht|wt|bmi).(.*)$")) %>%
  rename_with(year_to_age, matches("^(ht|wt|bmi).(.*)$")) %>%
  mutate(across(matches("^(ht|wt|bmi)_"), as.numeric),
         across(matches("^(ht|wt|bmi)_"), 
                ~ case_when(.x < 0 ~ NA_real_,
                            .x >= 7777 ~ NA_real_,
                            TRUE ~ .x)),
         ht_63 = ht_63*100) %>%
  # Survey Weight and Sex
  mutate(survey_weight = as.numeric(inf),
         female = sex - 1,) %>%
  # Social Class
  mutate(
    # Own Social Class
    own_class = factor(sep_levels[sc53u], sep_levels),
    own_class_ridit = make_ridit(own_class),
    # Own Education
    own_edu_level = names(educ_levels)[match(lhqr, educ_levels)] %>%
      factor(names(educ_levels)) %>% fct_rev(),
    own_edu_level_ridit = make_ridit(own_edu_level),
    # Father Class @ Age 4/11
    across(c(fsc50, fsc57),
           ~ as_factor(.x) %>%
             fct_recode(NULL = "Unknown") %>%
             as.numeric()),
    father_class = ifelse(!is.na(fsc50), fsc50, fsc57) %>%
      sep_levels[.] %>%
      factor(sep_levels),
    father_class_ridit = make_ridit(father_class),
    # Father's Education
    father_edu_years = clean_between(fagels, 10, 26) %>%
      truncate(10, 23) %>%
      subtract() %>%
      rescale(),
    father_edu_years_ridit = make_ridit(father_edu_years),
    father_edu_level = rec(fed, 
                           rec = "1 = 4 [Primary];
                     2 = 3 [Primary & FE or Higher];
                     3:5 = 2 [Secondary or Primary & FE or Higher];
                     6:8 = 1 [Secondary & FE or Higher]") %>%
      to_label(),
    father_edu_level_ridit = make_ridit(father_edu_level),
    # Mother's Education
    mother_edu_years = clean_between(magels, 10, 26) %>%
      truncate(10, 23) %>%
      subtract() %>%
      rescale(),
    mother_edu_years_ridit = make_ridit(mother_edu_years),
    mother_edu_level = rec(med, 
                           rec = "1 = 4 [Primary];
                     2 = 3 [Primary & FE or Higher];
                     3:5 = 2 [Secondary or Primary & FE or Higher];
                     6:8 = 1 [Secondary & FE or Higher]") %>%
      to_label(),
    mother_edu_level_ridit = make_ridit(mother_edu_level),
    # Condition of Child and Home
    condition_4 = factor(chilc50x, 0:2, c("Best", "Intermediate", "Worst")),
    condition_4_ridit = make_ridit(condition_4),
    housing_11 = ifelse(between(cohs57, 0, 5), 5 - cohs57, NA) %>%
      factor(),
    housing_11_ridit = make_ridit(housing_11),
  ) %>%
  # Choose Variables and Reshape
  select(id, survey_weight, female,
         matches("^(prs|bmi|wt|ht)_"), pc1:pc10,
         matches("^(own|mother|father)_")) %>%
  pivot_longer(matches("^(bmi|ht|wt)_"),
               names_to = c(".value", "age"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(age = as.numeric(age)) %>%
  arrange(id, age) %>%
  rename(height = ht, weight = wt) %>%
  mutate(bmi = ifelse(bmi > 80, NA, bmi),
         weight = ifelse(bmi >= 195, NA, weight)) %>%
  zap_formats() %>%
  zap_label() %>%
  zap_labels()


# Whole Body Measures
df_long <- df_raw %>%
  mutate(across(c(dxawbft09, dxawbln09), 
                ~ ifelse(.x %in% attr(.x, "labels"), NA, .x)),
         height = ifelse(htn09 >= 7777, NA, htn09),
         fat_ratio = ifelse(dxaangyr09 >= 7, NA, dxaangyr09),
         fat_mass = dxawbft09*0.001 / (height^2),
         lean_mass = dxawbln09*0.001 / (height^2), 
         across(fat_ratio:lean_mass, wtd_scale),
         age = 63) %>%
  select(id, fat_ratio, fat_mass, lean_mass, age) %>%
  full_join(df_long, by = c("id", "age"))

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
save(height_power, file = "Data/height_power.Rdata")

df_long <- df_long %>%
  left_join(height_power, by = "age") %>%
  mutate(bmi_corrected = weight / (height/100)^power) %>%
  select(-power) %>%
  filter(age != 6) %>%
  arrange(id, age) %>%
  relocate(id, age, survey_weight, prs_k, prs_v, prs_r_adult, prs_r_child)

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
  guides(color = "none")
ggsave("Images/height_power.png",
       width = 21, height = 21, units = "cm")


# Missing Social Class
df_raw %>%
  count(fsc50, fsc57) %>%
  mutate(across(1:2, ~ ifelse(between(.x, 10, 50), .x, NA))) %>%
  filter(is.na(fsc50)) %>%
  drop_na(fsc57) %>%
  count(wt = n)


# Correlation between SEP Measures
df_long %>%
  select(id, matches("ridit")) %>%
  distinct() %>%
  select(-id) %>%
  correlate() %>%
  shave(upper = TRUE) %>%
  stretch() %>%
  mutate(r = round(r, 2),
         across(x:y, fct_inorder)) %>%
  drop_na() %>%
  ggplot() +
  aes(x = x, y = y, label = r) +
  geom_text() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


df_long %>%
  select(id, matches("prs")) %>%
  distinct() %>%
  select(-id) %>%
  correlate() %>%
  shave(upper = TRUE)
