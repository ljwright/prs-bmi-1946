library(tidyverse)
library(haven)
library(glue)
library(broom)
library(gallimaufr)
library(magrittr)
library(quantreg)

rm(list = ls())

# 1. Load Data ----
df <- read_dta("C:/Users/liamj/OneDrive - University College London/prs_bmi/analysis/output/46c_bmi_ht_bmiprs.dta")  %>%
  mutate(id = row_number(), .before = 1)

sexes <- list(female = "women",
              male = "men",
              all = c("women", "men"))


# 2. Linear Regression Models ----
id_all <- df %>%
  select(id, matches("logbmi")) %>%
  drop_na() %>%
  select(id)

df_bmi <- df %>%
  select(id, matches("^(logbmi|bmi)\\d+$"), -bmi09) %>%
  rename_with(~ str_replace(.x, "bmi", "bmi_")) %>% # names()
  pivot_longer(-id, names_to = c(".value", "age"),
               names_pattern = "(.*)_(.*)") %>%
  rename(bmi_ln = logbmi) %>%
  mutate(bmi = as.numeric(bmi),
         age = as.numeric(age)) %>%
  group_by(age) %>%
  mutate(bmi_std = wtd_scale(bmi)) %>%
  ungroup()

df_prs <- df %>%
  select(id, sex, matches("^zprs_")) %>%
  mutate(sex = as_factor(sex)) %>%
  drop_na()

mod_specs <- expand_grid(age = unique(df_bmi$age),
                         sex = names(sexes),
                         prs = str_subset(names(df_prs), "zprs"),
                         bmi_var = str_subset(names(df_bmi), "bmi"),
                         obs = c("cc", "obs")) %>%
  mutate(spec_id = row_number(), .before = 1)

get_res <- function(spec_id){
  spec <- mod_specs %>%
    filter(spec_id == !!spec_id) %>%
    as.list()
  
  # Select Data
  df_mod <- df_bmi %>%
    left_join(df_prs, by = "id") %>%
    filter(sex %in% sexes[[!!spec$sex]],
           age == !!spec$age) %>%
    rename(prs = all_of(spec$prs))
  
  if (spec$obs == "cc"){
    df_mod <- inner_join(df_mod, id_all, by = "id")
  }
  
  # Run Model
  mod_form <- glue("{spec$bmi_var} ~ prs")
  if (spec$sex == "all") mod_form <- glue("{mod_form} + sex")
  
  mod <- as.formula(mod_form) %>%
    lm(df_mod)
  
  tidy(mod, conf.int = TRUE) %>%
    filter(term == "prs") %>%
    select(beta = estimate, p = p.value, lci = conf.low, uci = conf.high) %>%
    bind_cols(glance(mod) %>% select(r2 = 1, n = 12))
}

gwas_dict <- c(zprs_k = "Khera et al. (2019)", 
               zprs_r = "Richardson et al. (2020)", 
               zprs_v = "Vogelezang et al. (2020)")

df_res <- mod_specs %>%
  mutate(map_dfr(spec_id, get_res)) %>%
  select(-spec_id)

plot_res <- function(obs, bmi_var, save_p = FALSE){
  p <- df_res %>%
    filter(obs == !!obs, bmi_var == !!bmi_var) %>%
    pivot_longer(c(beta, r2), names_to = "stat", values_to = "beta") %>%
    mutate(across(c(lci, uci), ~ ifelse(stat == "r2", NA, .x)),
           age = factor(age) %>% fct_rev(),
           sex = str_to_title(sex),
           stat = ifelse(stat == "beta", "Effect Size", "R-Squared"),
           gwas_clean = factor(gwas_dict[prs], gwas_dict)) %>%
    ggplot() +
    aes(x = age, y = beta, ymin = lci, ymax = uci, color = sex) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(gwas_clean ~ stat, scales = "free_x", switch = "y") +
    coord_flip() +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom") +
    labs(x = NULL, y = NULL, color = NULL)
  
  if (save_p){
    glue("Images/res_{obs}_{bmi_var}.png") %>%
      ggsave(p, height = 29.7, width = 21, units = "cm")
  }
  
  return(p)
}

df_res %>%
  distinct(obs, bmi_var) %$%
  map2(obs, bmi_var, plot_res, TRUE)

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
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0)) +
    labs(x = "Polygenic Risk Score", y = "BMI", color = "Age") +
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
  theme_minimal() +
  labs(x = "BMI", y = "Density")
ggsave("Images/density.png",
       height = 21, width = 29.7, units = "cm")


# 5. Quantile Regression ----
get_quant <- function(spec_id){
  spec <- mod_specs %>%
    filter(spec_id == !!spec_id) %>%
    as.list()
  
  # Select Data
  df_mod <- df_bmi %>%
    left_join(df_prs, by = "id") %>%
    filter(sex %in% sexes[[!!spec$sex]],
           age == !!spec$age) %>%
    rename(prs = all_of(spec$prs))
  
  if (spec$obs == "cc"){
    df_mod <- inner_join(df_mod, id_all, by = "id")
  }
  
  # Run Model
  mod_form <- glue("{spec$bmi_var} ~ prs")
  if (spec$sex == "all") mod_form <- glue("{mod_form} + sex")
  
  get_rq <- function(tau){
    as.formula(mod_form) %>%
      rq(tau = tau, data = df_mod) %>%
      tidy(conf.int = TRUE) %>%
      filter(term == "prs") %>%
      select(beta = estimate, lci = conf.low, uci = conf.high)
  }
  
  tibble(tau = 1:9/10) %>%
    mutate(map_dfr(tau, get_rq))
}

df_quant <- mod_specs %>%
  mutate(res = map(spec_id, get_quant)) %>%
  unnest(res) %>%
  select(-spec_id)

plot_quant <- function(obs, bmi_var, save_p = FALSE){
  p <- df_quant %>%
    filter(obs == !!obs, bmi_var == !!bmi_var) %>%
    mutate(age = factor(age) %>% ordered(),
           sex = str_to_title(sex),
           gwas_clean = factor(gwas_dict[prs], gwas_dict))%>%
    ggplot() +
    aes(x = tau, y = beta, ymin = lci, ymax = uci,
        color = age, fill = age) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(gwas_clean ~ sex, scales = "free", switch = "y") +
    geom_ribbon(color = NA, alpha = 0.1) +
    geom_line() +
    scale_x_continuous(breaks = 1:9/10, labels = glue("{1:9*10}th")) +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Percentile", y = NULL, color = "Age", fill = "Age")
  
  if (save_p){
    glue("Images/quant_{obs}_{bmi_var}.png") %>%
      ggsave(p, height = 29.7, width = 21, units = "cm")
  }
  
  return(p)
}

df_quant %>%
  distinct(obs, bmi_var) %$%
  map2(obs, bmi_var, plot_quant, TRUE)
