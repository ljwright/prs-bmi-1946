library(tidyverse)
library(haven)
library(glue)
library(broom)
library(gallimaufr)
library(magrittr)
library(quantreg)
library(splines)
library(furrr)
library(tictoc)
library(margins)

rm(list = ls())

plan(multisession, workers = 8)

# 1. Load Data ----
load("Data/df_long.Rdata")

df_reg <- df_long %>%
  group_by(age) %>%
  mutate(across(c(bmi, bmi_corrected),
                list(std = wtd_scale,
                     ln = log,
                     rank = ~ percent_rank(.x)*100))) %>%
  ungroup() %>%
  rename(bmi_raw = bmi, bmi_corrected_raw = bmi_corrected) %>%
  drop_na(matches("^prs"))

rm(df_long)

# Model Objects
ages <- unique(df_reg$age)

sexes <- list(female = 1,
              male = 0,
              all = c(0, 1))

pcs <- paste0("pc", 1:10)

dep_vars <- names(df_reg) %>%
  str_subset("^(bmi_|height|weight|fat_|lean)")

sep_vars <- str_subset(names(df_reg), "(edu_|class|ridit)")
all_sep <- sep_vars %>% 
  str_subset("(edu_level|class).*ridit", TRUE) %>% 
  str_subset("edu_years$", TRUE)

alt_dep <- c("fat_ratio", "fat_mass", "lean_mass")

mod_specs <- expand_grid(age = ages,
                         sex = names(sexes),
                         prs_var = str_subset(names(df_reg), "prs"),
                         dep_var = dep_vars) %>%
  filter(!(dep_var %in% alt_dep) |
           (dep_var %in% alt_dep & age == 63)) %>%
  mutate(spec_id = row_number(), .before = 1)

sep_specs <- mod_specs %>%
  expand_grid(sep_var = c(as.list(sep_vars),
                          list(all_sep)))


# 2. Model Functions ----
get_spec <- function(spec_id){
  spec <- mod_specs %>%
    filter(spec_id == !!spec_id) %>%
    as.list()
  
  return(spec)
}

get_df <- function(spec_id){
  spec <- get_spec(spec_id)
  
  df_reg %>%
    rename(prs = all_of(spec$prs_var), 
           dep_var = all_of(spec$dep_var)) %>%
    filter(female %in% sexes[[!!spec$sex]],
           age == !!spec$age)
}

get_form <- function(spec_id, covars = c("prs", pcs)){
  spec <- get_spec(spec_id)
  
  mod_form <- glue_collapse(covars, ' + ') %>%
    paste("dep_var ~", .)
  
  if (spec$sex == "all") mod_form <- glue("{mod_form} + female")
  
  return(mod_form)
}

get_ci <- function(estimates){
  quantile(estimates, c(.5, .025, .975)) %>%
    as_tibble_row() %>%
    rename(beta = 1, lci = 2, uci = 3)
}

get_furrr <- function(df_specs, func){
  df_specs %>%
    sample_frac() %>%
    mutate(res = future_map(spec_id, func, 
                            .progress = TRUE,
                            .options = furrr_options(seed = NULL))) %>%
    unnest(res) %>%
    arrange(spec_id)
}

get_furrr2 <- function(df_specs, var, func){
  df_specs %>%
    sample_frac() %>%
    mutate(res = future_map2(spec_id, {{ var }}, func, 
                             .progress = TRUE,
                             .options = furrr_options(seed = NULL))) %>%
    unnest(res) %>%
    arrange(spec_id)
}


# 3. Main Regressions ----
get_lm <- function(spec_id){
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id) %>%
    select(id, dep_var, prs, female, pc1:pc10) %>%
    drop_na()
  
  # Main Effect
  mod <- get_form(spec_id) %>%
    as.formula() %>%
    lm(df_mod)
  
  mod_coef <- broom::tidy(mod, conf.int = TRUE) %>%
    filter(term == "prs") %>%
    select(term, beta = estimate, lci = conf.low, uci = conf.high)
  
  # Incremental R2
  get_r2 <- function(df_m){
    r2_prs <- get_form(spec_id) %>%
      as.formula() %>%
      lm(df_m) %>%
      broom::glance() %>%
      pull(r.squared)
    
    r2_other <- get_form(spec_id, pcs) %>%
      as.formula() %>%
      lm(df_m) %>%
      broom::glance() %>%
      pull(r.squared)
    
    r2_prs - r2_other
  }
  
  get_boot <- function(boot){
    set.seed(boot)
    
    sample_frac(df_mod, replace = TRUE) %>%
      get_r2()
  }
  
  map_dbl(1:500, get_boot) %>%
    get_ci() %>%
    mutate(term = "r2",
           beta = get_r2(df_mod),
           .before = 1) %>%
    bind_rows(mod_coef)
}

tic()
res_prs <- get_furrr(mod_specs, get_lm)
toc()
save(res_prs, file = "Data/res_prs.Rdata")
load("Data/res_prs.Rdata")


# 4. Attrition Regressions ----
attrit_sample <- function(spec_id, age_from){
  spec <- get_spec(spec_id)
  
  df_reg %>%
    filter(age >= age_from) %>%
    drop_na(all_of(c(spec$prs_var, spec$dep_var))) %>%
    count(id) %>%
    filter(n == length(ages[ages >= age_from])) %>%
    pull(id)
}

#  Main Effect
get_attrit <- function(spec_id, age_from){
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id)
  
  if (age_from > 0){
    df_mod <- df_mod %>%
      filter(id %in% attrit_sample(!!spec_id, !!age_from))
  }
  
  get_form(spec_id) %>%
    as.formula() %>%
    lm(df_mod) %>%
    broom::tidy(conf.int = TRUE) %>%
    filter(term == "prs") %>%
    select(beta = 2, p = 5, lci = 6, uci = 7)
}

tic()
res_attrit <- expand_grid(mod_specs, age_from = c(0, ages)) %>%
  filter(age >= age_from,
         str_detect(dep_var, "bmi")) %>%
  get_furrr2(age_from, get_attrit)
toc()

# R-Squared
get_attrit_r2 <- function(spec_id, age_from){
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id)
  
  if (age_from > 0){
    df_mod <- df_mod %>%
      filter(id %in% attrit_sample(!!spec_id, !!age_from))
  }
  
  get_r2 <- function(df_m){
    r2_prs <- get_form(spec_id) %>%
      as.formula() %>%
      lm(df_m) %>%
      broom::glance() %>%
      pull(r.squared)
    
    r2_other <- get_form(spec_id, pcs) %>%
      as.formula() %>%
      lm(df_m) %>%
      broom::glance() %>%
      pull(r.squared)
    
    r2_prs - r2_other
  }
  
  get_boot <- function(boot){
    set.seed(boot)
    
    sample_frac(df_mod, replace = TRUE) %>%
      get_r2()
  }
  
  map_dbl(1:500, get_boot) %>%
    get_ci() %>%
    mutate(beta = get_r2(df_mod))
}

tic()
res_attrit_r2 <- expand_grid(mod_specs, age_from = c(0, ages)) %>%
  filter(age >= age_from,
         dep_var == "bmi_raw",
         sex == "all") %>%
  get_furrr2(age_from, get_attrit_r2)
toc()
save(res_attrit, res_attrit_r2, file = "Data/res_attrit.Rdata")
load("Data/res_attrit.Rdata")


# 5. Quantile Regression ----
get_quant <- function(spec_id){
  df_mod <- get_df(spec_id)
  
  mod_form <- get_form(spec_id)
  
  as.formula(mod_form) %>%
    rq(tau = 1:9/10, data = df_mod) %>%
    broom::tidy(conf.int = TRUE) %>%
    filter(term == "prs") %>%
    select(tau, beta = estimate, lci = conf.low, uci = conf.high)
  
}

tic()
res_quant <- mod_specs %>%
  filter(str_detect(dep_var, "(raw|std|ln)")) %>%
  get_furrr(get_quant)
toc()

# 5. SEP Additive ----
get_sep <- function(spec_id, sep_vars){
  
  df_mod <- get_df(spec_id) %>%
    select(id, dep_var, prs, all_of(sep_vars), female, pc1:pc10) %>%
    drop_na()
  
  reg_ex <- sep_vars %>%
    glue_collapse("|") %>%
    paste0("^(", ., ")")
  
  run_sep <- function(covars){
    mod <- get_form(spec_id, covars) %>%
      as.formula() %>%
      lm(df_mod)
    
    broom::tidy(mod, conf.int = TRUE) %>%
      filter(str_detect(term, reg_ex)) %>%
      select(term, beta = 2, p = 5, lci = 6, uci = 7) %>%
      bind_cols(broom::glance(mod) %>%
                  select(r2 = r.squared, n = nobs))
  }
  
  get_r2 <- function(covars){
    get_form(spec_id, covars) %>%
      as.formula() %>%
      lm(df_mod) %>%
      broom::glance() %>%
      pull(r.squared)
  }
  
  r_prs <- get_r2(c("prs", pcs))
  r_sex <- get_r2("1")
  
  bind_rows(bivar = run_sep(sep_vars),
            adjust = run_sep(c("prs", pcs, sep_vars)),
            .id = "mod") %>%
    mutate(r_prs = r_prs, r_sex = r_sex,
           r_prs_diff = r2 - r_prs,
           r_sex_diff = r2 - r_sex)
  
}

tic()
res_sep <- sep_specs %>%
  get_furrr2(sep_var, get_sep)
toc()


# 7. SEP Multiplicative ----
get_mult <- function(spec_id, sep_var){
  
  sep_vals <- df_reg %>%
    rename(sep_var = all_of(!!sep_var)) %>%
    count(sep_var) %>%
    drop_na() %>%
    pull(1) %>%
    c(0, ., 1) %>%
    unique()
  
  df_mod <- get_df(spec_id) %>%
    rename(sep_var = all_of(!!sep_var))
  
  mod <- get_form(spec_id, c("prs*sep_var", pcs)) %>%
    as.formula() %>%
    lm(df_mod)
  
  mrg <- margins(mod, at = list(sep_var = sep_vals), variables = "prs") %>%
    broom::tidy(conf.int = TRUE) %>%
    select(sep_values = 3, beta = 4, p = 7, lci = 8, uci = 9) %>%
    distinct() %>%
    list()
  
  res <- broom::tidy(mod, conf.int = TRUE) %>%
    filter(term == "prs:sep_var") %>%
    select(beta = 2, p = 5, lci = 6, uci = 7) %>%
    bind_cols(glance(mod) %>% select(r2 = r.squared, n = nobs)) %>%
    list()
  
  tibble(margins = mrg, mod = res)
  
}

tic()
res_mult <- sep_specs %>%
  filter(map_dbl(sep_var, length) == 1) %>%
  unnest(sep_var) %>%
  filter(sep_var %in% str_subset(sep_vars, "(ridit|edu_years)")) %>%
  get_furrr2(sep_var, get_mult)
toc()


# 8. Save Objects ----
save(res_prs, res_attrit, res_attrit_r2,
     res_quant, res_sep, res_mult, 
     mod_specs, get_ci,
     file = "Data/regression_results.Rdata")
