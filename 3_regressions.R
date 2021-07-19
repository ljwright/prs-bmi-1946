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

plan(multisession, workers = 4)

# 1. Load Data ----
load("Data/df_long.Rdata")

df_reg <- df_long %>%
  group_by(age) %>%
  mutate(across(c(bmi, bmi_corrected),
                list(std = wtd_scale,
                     ln = log,
                     rank = ~ percent_rank(.x)*100))) %>%
  ungroup() %>%
  rename(bmi_raw = bmi, bmi_corrected_raw = bmi_corrected)

rm(df_long)

# Model Objects
ages <- unique(df_reg$age)

sexes <- list(female = 1,
              male = 0,
              all = c(0, 1))

dep_vars <- names(df_reg) %>%
  str_subset("^(bmi|height|weight|fat|lean)")

alt_dep <- c("fat_ratio", "fat_mass", "lean_mass")

mod_specs <- expand_grid(age = ages,
                         sex = names(sexes),
                         prs_var = str_subset(names(df_reg), "prs"),
                         dep_var = dep_vars) %>%
  filter(!(dep_var %in% alt_dep) |
           (dep_var %in% alt_dep & age == 63)) %>%
  mutate(spec_id = row_number(), .before = 1)

sep_specs <- mod_specs %>%
  expand_grid(sep_var = c("medu", "medu_ridit", "sep_ridit"))


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

get_form <- function(spec_id, covars = "prs"){
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
    mutate(res = future_map(spec_id, func, 
                            .progress = TRUE,
                            .options = furrr_options(seed = NULL))) %>%
    unnest(res)
}

get_furrr2 <- function(df_specs, var, func){
  df_specs %>%
    mutate(res = future_map2(spec_id, {{ var }}, func, 
                             .progress = TRUE,
                             .options = furrr_options(seed = NULL))) %>%
    unnest(res)
}


# 3. Main Regressions ----
get_lm <- function(spec_id){
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id) %>%
    select(id, dep_var, prs, female) %>%
    drop_na()
  
  mod_form <- get_form(spec_id)
  
  get_boot <- function(boot){
    set.seed(boot)
    
    df_m <- sample_frac(df_mod, replace = TRUE)
    
    mod <- as.formula(mod_form) %>%
      lm(df_m)
    
    r2 <- broom::glance(mod)[[1]]
    if (spec$sex == "all"){
      r2_base <- get_form(spec_id, "1") %>%
        lm(df_m) %>%
        broom::glance() %>%
        pull(1)
      r2 <- r2 - r2_base
    }
    
    c(coef(mod)["prs"], r2 = r2) %>%
      enframe(name = "term", value = "estimate")
  }
  
  map_dfr(1:500, get_boot, .id = "boot") %>%
    group_by(term) %>%
    summarise(get_ci(estimate)) %>%
    mutate(n = nrow(df_mod))
}

tic()
res_prs <- get_furrr(mod_specs, get_lm)
toc()


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
get_sep <- function(spec_id, sep_var){
  
  df_mod <- get_df(spec_id) %>%
    rename(sep_var = all_of(!!sep_var)) %>%
    select(id, dep_var, prs, sep_var, female) %>%
    drop_na()
  
  run_sep <- function(covars){
    mod <- get_form(spec_id, covars) %>%
      as.formula() %>%
      lm(df_mod)
    
    broom::tidy(mod, conf.int = TRUE) %>%
      filter(term == "sep_var") %>%
      select(beta = 2, p = 5, lci = 6, uci = 7) %>%
      bind_cols(broom::glance(mod) %>%
                  select(r2 = 1, n = 12))
  }
  
  r_comp <- get_form(spec_id) %>%
    as.formula() %>%
    lm(df_mod) %>%
    broom::glance(mod) %>%
    pull(1)
  
  bind_rows(bivar = run_sep("sep_var"),
            adjust = run_sep(c("prs", "sep_var")),
            .id = "mod") %>%
    mutate(r_comp = r_comp,
           r_diff = r2 - r_comp)
  
}

tic()
res_sep <- sep_specs %>%
  get_furrr2(sep_var, get_sep)
toc()

# Bootstraps for R2
# get_sep_boot <- function(spec_id, sep_var){
#   
#   spec <- get_spec(spec_id)
#   
#   df_mod <- get_df(spec_id) %>%
#     rename(sep_var = all_of(!!sep_var)) %>%
#     select(id, dep_var, prs, sep_var, female) %>%
#     drop_na()
#   
#   get_boot <- function(boot){
#     set.seed(boot)
#     
#     df_m <- sample_frac(df_mod, replace = TRUE)
#     
#     run_sep <- function(covars){
#       mod <- get_form(spec_id, covars) %>%
#         as.formula() %>%
#         lm(df_m)
#       
#       coefs <- coef(mod)
#       
#       c(coefs["sep_var"],
#         r2 = broom::glance(mod)[[1]]) %>%
#         enframe(name = "term", value = "estimate")
#     }
#     
#     res <- bind_rows(bivar = run_sep("sep_var"),
#                      adjust = run_sep(c("prs", "sep_var")),
#                      .id = "mod")
#     
#     r2_prs <- get_form(spec_id, "prs") %>%
#       as.formula() %>%
#       lm(df_m) %>%
#       broom::glance() %>%
#       pull(1)
#     
#     res %>%
#       uncount(ifelse(term == "r2", 2, 1), .id = "id") %>%
#       mutate(term = ifelse(id == 2, "r2_diff", term),
#              estimate = ifelse(id == 2, estimate - r2_prs, estimate))
#   }
#   
#   map_dfr(1:500, get_boot, .id = "boot") %>%
#     group_by(mod, term) %>%
#     summarise(get_ci(estimate),
#               .groups = "drop") %>%
#     mutate(n = nrow(df_mod))
#   
# }
# 
# tic()
# res_sep_boot <- sep_specs %>%
#   filter(sep_var == "sep_ridit") %>%
#   get_furrr2(sep_var, get_sep_boot)
# toc()

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
  
  mod <- get_form(spec_id, "prs*sep_var") %>%
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
    bind_cols(glance(mod) %>% select(r2 = 1, n = 12)) %>%
    list()
  
  tibble(margins = mrg, mod = res)
  
}

tic()
res_mult <- get_furrr2(sep_specs, sep_var, get_mult)
toc()

# 8. SEP Quantile Regression ----
get_sep_quant <- function(spec_id, sep_var){
  
  df_mod <- get_df(spec_id) %>%
    rename(sep_var = all_of(!!sep_var))
  
  mod_form <- get_form(spec_id, "prs*sep_var")
  
  as.formula(mod_form) %>%
    rq(tau = 1:9/10, data = df_mod) %>%
    broom::tidy(conf.int = TRUE) %>%
    filter(term == "prs:sep_var") %>%
    select(tau, beta = estimate, lci = conf.low, uci = conf.high)
  
}

tic()
res_sep_quant <- sep_specs %>%
  filter(str_detect(dep_var, "(raw|std|ln)")) %>%
  get_furrr2(sep_var, get_sep_quant)
toc()


# 10. Save Objects ----
save(res_prs, res_attrit, res_quant,
     res_sep, res_mult, res_sep_quant,
     mod_specs, get_ci,
     file = "Data/regression_results.Rdata")
