library(tidyverse)
library(haven)
library(glue)
library(broom)
library(gallimaufr)
library(magrittr)

rm(list = ls())

# KERNEL DENSITY FOR EACH SAMPLE AND SEX
# DENSITY OF PRS BY FOLLOW-UP
# DENSITY OF BMI BY WHETHER PRS OBSERVED

# 1. Load Data ----
load("Data/df_long.Rdata")


# 2. Attrition ----
id_follow <- unique(df_long$age) %>%
  set_names(., .) %>%
  map(~ filter(df_long, age == .x) %>%
        drop_na(bmi) %>%
        distinct(id) %>%
        pull(id)) %>%
  c(list(obs = unique(df_long$id),
         cc = df_long %>%
           select(id, age, bmi) %>%
           pivot_wider(names_from = age, values_from = bmi) %>%
           drop_na() %>% 
           pull(id)))


df_long %>%
  mutate(miss_prs = if_any(matches("prs"), is.na) %>%
           if_else("Observed PRS", "Missing PRS") %>%
           factor()) %>%
  ggplot() +
  aes(x = bmi, color = miss_prs, fill = miss_prs) +
  facet_wrap(~ age, scales = "free") +
  geom_density(alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "BMI", y = "Density", color = NULL, fill = NULL)

df_long %>%
  drop_na(bmi) %>%
  mutate(age = factor(age)) %>%
  ggplot() +
  aes(x = prs_k) +
  facet_wrap(~ age, scales = "free") +
  geom_density(data = distinct(df_long, id, prs_k),
               color = "grey70", fill = "grey70", alpha = 0.4) +
  geom_density(aes(color = age, fill = age), alpha = 0.4) +
  guides(color = FALSE, fill = FALSE) +
  theme_bw() +
  labs(x = "Polygenic Risk Score", y = "Density")


get_lm <- function(obs_age, reg_age){
  df_long %>%
    filter(id %in% id_follow[[!!obs_age]]) %>%
    drop_na(prs_k, bmi) %>%
    distinct(id) %>%
    left_join(df_long, by = "id") %>%
    filter(age == !!reg_age) %>%
    lm(bmi ~ prs_k + female, .) %>%
    tidy(conf.int = TRUE) %>%
    filter(term == "prs_k") %>%
    select(beta = 2, lci = 6, uci = 7)
}

attrit_lm <- expand_grid(obs_age = names(id_follow),
                         reg_age = unique(df_long$age)) %>%
  mutate(map2_dfr(obs_age, reg_age, get_lm),
         across(matches("age"), ordered),
         obs_age = factor(obs_age, names(id_follow)) %>%
           fct_recode("Observed" = "obs", "Complete Cases" = "cc"))

ggplot(attrit_lm) +
  aes(x = reg_age, y = beta, ymin = lci, ymax = uci) +
  facet_wrap(~ obs_age) +
  geom_hline(yintercept = 0) +
  geom_line(data = rename(attrit_lm, obs = obs_age),
            aes(group = obs), color = "grey70") +
  geom_ribbon(aes(fill = obs_age, group = obs_age), color = NA, alpha = 0.2) +
  geom_line(aes(color = obs_age, group = obs_age), size = 1) +
  theme_minimal() +
  labs(x = "Age", y = "Marginal Effect") +
  guides(color = FALSE, fill = FALSE)


# 2. Model Objects ----
id_all <- df %>%
  select(id, matches("logbmi")) %>%
  drop_na() %>%
  select(id)

df_bmi <- df %>%
  select(id, matches("^(logbmi|bmi)\\d+$"), -bmi09) %>%
  rename_with(~ str_replace(.x, "bmi", "bmi_")) %>%
  pivot_longer(-id, names_to = c(".value", "age"),
               names_pattern = "(.*)_(.*)") %>%
  rename(bmi_ln = logbmi) %>%
  mutate(bmi = as.numeric(bmi),
         age = as.numeric(age)) %>%
  group_by(age) %>%
  mutate(bmi_std = wtd_scale(bmi)) %>%
  ungroup()

df_prs <- df %>%
  select(id, sex, matches("^zprs_"), fsc4_ridit) %>%
  as_factor() %>%
  drop_na(matches("zprs"))

sexes <- list(female = "women",
              male = "men",
              all = c("women", "men"))

mod_specs <- expand_grid(age = unique(df_bmi$age),
                         sex = names(sexes),
                         prs = str_subset(names(df_prs), "zprs"),
                         bmi_var = str_subset(names(df_bmi), "bmi"),
                         obs = c("cc", "obs")) %>%
  mutate(spec_id = row_number(), .before = 1)

save(df, id_all, df_bmi, df_prs, sexes, mod_specs,
     file = "Data/model_objects.Rdata")

# 3. Scatter Plots ----
plot_scatter <- function(prs, save_p = FALSE){
  df_scat <- df_bmi %>%
    left_join(df_prs, by = "id") %>%
    select(id, age, bmi, matches("zprs")) %>%
    pivot_longer(matches("zprs"), names_to = "prs", values_to = "prs_value") %>%
    drop_na() %>%
    inner_join(id_all, by = "id") %>%
    mutate(age = factor(age) %>% ordered(),
           bmi = as.numeric(bmi)) %>%
    filter(prs == !!prs)
  
  p <- ggplot(df_scat) +
    aes(x = bmi, y = prs_value) +
    facet_wrap(~ age, ncol = 4) +
    coord_flip() +
    geom_jitter(data = select(df_scat, -age), color = "grey70", alpha = 0.2) +
    geom_jitter(aes(color = age), alpha = 0.2) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0)) +
    labs(x = "BMI", y = "Polygenic Risk Score", color = "Age") +
    guides(color = FALSE)
  
  if (save_p){
    glue("Images/scatter_{prs}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

mod_specs %>%
  distinct(prs) %$%
  walk(prs, plot_scatter, TRUE)


# 4. Kernel Density ----
df_dens <- df_bmi %>%
  select(id, age, bmi) %>%
  drop_na() %>%
  inner_join(id_all, by = "id") %>%
  mutate(age = factor(age) %>% ordered(),
         age_f = age,
         bmi = as.numeric(bmi))

ggplot(df_dens) +
  aes(x = bmi) +
  facet_wrap(~ age, ncol = 4) +
  geom_density(data = select(df_dens, -age), aes(group = age_f),
               color = "grey70", fill = "grey70", alpha = 0.4) +
  geom_density(aes(color = age, fill = age), alpha = 0.7) +
  guides(color = FALSE, fill = FALSE) +
  theme_bw() +
  labs(x = "BMI", y = "Density")
ggsave("Images/density.png",
       height = 21, width = 29.7, units = "cm")

# 4. PRS x SEP
gwas_dict <- c(zprs_k = "Khera et al. (2019)", 
               zprs_r = "Richardson et al. (2020)", 
               zprs_v = "Vogelezang et al. (2020)")

df_fsc4 <- df %>%
  distinct(fsc4, fsc4_ridit) %>%
  mutate(fsc4 = case_when(fsc4 == 0 ~ "I",
                          fsc4 == 1 ~ "II",
                          fsc4 == 2 ~ "III",
                          fsc4 == 3 ~ "IV",
                          fsc4 == 4 ~ "V Skilled",
                          fsc4 == 5 ~ "V Unskilled",
                          TRUE ~ NA_character_) %>%
           fct_reorder(fsc4) %>% 
           factor(ordered = TRUE)) %>%
  drop_na() %>%
  left_join(df_prs, by = "fsc4_ridit") %>%
  pivot_longer(matches("zprs"), names_to = "prs_var", values_to = "prs") %>%
  mutate(fsc4_f = fsc4,
         gwas_clean = factor(gwas_dict[prs_var], gwas_dict)) %>%
  select(id, fsc4, fsc4_f, gwas_clean, prs)

ggplot(df_fsc4) +
  aes(x = prs, group = fsc4_f) +
  facet_grid(gwas_clean ~ fsc4, switch = "y") +
  geom_density(data = select(df_fsc4, -fsc4),
               color = "grey70", fill = "grey70", alpha = 0.4) +
  geom_density(aes(color = fsc4, fill = fsc4), alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0)) +
  labs(x = "Polygenic Risk Score", y = NULL,
       color = "Social Class", fill = "Social Class")
ggsave("Images/density_sep.png",
       height = 21, width = 29.7, units = "cm")

df_fsc4 %>%
  mutate(fsc4 = factor(fsc4, ordered = FALSE)) %>%
  nest(data = -gwas_clean) %>%
  mutate(res = map(data,
                   ~ lm(prs ~ fsc4, .x) %>%
                     tidy(conf.int = TRUE))) %>%
  unnest(res) %>%
  select(-data) %>%
  filter(str_detect(term, "fsc4")) %>%
  mutate(term = str_replace(term, "fsc4", "") %>%
           ordered(levels(df_fsc4$fsc4)) %>%
           fct_rev()) %>%
  ggplot() +
  aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high) +
  facet_wrap(~ gwas_clean) +
  geom_hline(yintercept = 0) +
  geom_pointrange() +
  theme_minimal() +
  coord_flip() +
  labs(x = NULL, y = "Difference in Polygenic Risk Score")
ggsave("Images/diff_prs_sep.png",
       height = 16, width = 21, units = "cm")

# 5. Probability of Superiority ----
get_superior <- function(left, right, gwas_clean){
  set.seed(1)
  
  get_prs <- function(fsc4, gwas_clean){
    df_fsc4 %>%
      filter(fsc4 == !!fsc4,
             gwas_clean == !!gwas_clean) %>%
      sample_n(1000, TRUE) %>%
      pull(prs)
  }
  
  prs_left <- get_prs(left, gwas_clean)
  prs_right <- get_prs(right, gwas_clean)
  
  tibble(diff = prs_left - prs_right)
}


df_prob <- df_fsc4 %$%
  expand_grid(left = unique(fsc4),
              right = unique(fsc4), 
              gwas_clean = unique(gwas_clean)) %>%
  mutate(res = pmap(list(left, right, gwas_clean), get_superior)) %>%
  unnest(res) %>%
  filter(left != right)

df_prob_x <- df_prob %>%
  group_by(across(-diff)) %>%
  summarise(p = sum(diff > 0)/n(),
            Mean = mean(diff),
            Median = median(diff),
            .groups = "drop") %>%
  pivot_longer(c(p, Mean, Median)) %>%
  filter(gwas_clean == "Khera et al. (2019)") %>%
  mutate(x = Inf, y = Inf)

df_prob %>%
  filter(gwas_clean == "Khera et al. (2019)") %>%
  ggplot() +
  aes(x = diff) +
  facet_grid(left ~ right) +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_vline(data = filter(df_prob_x, name != "p"),
             aes(xintercept = value, linetype = name,
                 color = name)) +
  geom_label(data = filter(df_prob_x, name == "p"),
             aes(x = x, y = y, label = value)) +
  geom_density(color = "grey60", fill = "grey60", alpha = 0.4) +
  theme_minimal() +
  labs(x = "Differences in Polygenic Risk Scores",
       y = "Density")
