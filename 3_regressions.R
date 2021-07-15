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

df_long <- df_long %>%
  group_by(age) %>%
  mutate(across(c(bmi, bmi_corrected),
                list(std = wtd_scale,
                     ln = log,
                     rank = ~ percent_rank(.x)*100))) %>%
  ungroup() %>%
  rename(bmi_raw = bmi, bmi_corrected_raw = bmi_corrected) %>%
  select(-height, -weight)

ages <- unique(df_long$age)

sexes <- list(female = 1,
              male = 0,
              all = c(0, 1))

mod_specs <- expand_grid(age = ages,
                         sex = names(sexes),
                         prs_var = str_subset(names(df_long), "prs"),
                         dep_var = str_subset(names(df_long), "^bmi")) %>%
  mutate(spec_id = row_number(), .before = 1)


# 2. Model Functions ----
get_spec <- function(spec_id){
  spec <- mod_specs %>%
    filter(spec_id == !!spec_id) %>%
    as.list()
  
  return(spec)
}

get_df <- function(spec_id){
  spec <- get_spec(spec_id)
  
  df_long %>%
    rename(prs = all_of(spec$prs_var), 
           dep_var = all_of(spec$dep_var)) %>%
    filter(female %in% sexes[[!!spec$sex]],
           age == !!spec$age) %>%
    select(id, dep_var, prs, sep_ridit, female)
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


# 3. Main Regressions ----
get_lm <- function(spec_id){
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id) %>%
    select(-sep_ridit) %>%
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
attrit_sample <- function(age_from, spec_id){
  spec <- get_spec(spec_id)
  
  df_long %>%
    filter(age >= age_from) %>%
    drop_na(all_of(c(spec$prs_var, spec$dep_var))) %>%
    count(id) %>%
    filter(n == length(ages[ages >= age_from])) %>%
    pull(id)
}

get_attrit <- function(age_from, spec_id){
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id)
  
  if (age_from > 0){
    df_mod <- df_mod %>%
      filter(id %in% attrit_sample(!!age_from, !!spec_id))
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
  filter(age >= age_from) %>%
  mutate(future_map2_dfr(age_from, spec_id, get_attrit, 
                         .progress = TRUE,
                         .options = furrr_options(seed = NULL)),
         age_from = factor(age_from) %>%
           fct_recode("Observed" = "0", "Complete Cases" = "2"))
toc()


# 5. Quantile Regression ----
get_quant <- function(spec_id){
  df_mod <- get_df(spec_id)
  
  mod_form <- get_form(spec_id)
  
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

tic()
res_quant <- mod_specs %>%
  filter(str_detect(dep_var, "(raw|std|ln)")) %>%
  get_furrr(get_quant)
toc()

# 5. SEP Additive ----
get_sep <- function(spec_id){
  
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id) %>% 
    drop_na()
  
  mod_form <- get_form(spec_id)
  
  get_boot <- function(boot){
    set.seed(boot)
    
    df_m <- sample_frac(df_mod, replace = TRUE)
    
    mod <- as.formula(mod_form) %>%
      lm(df_m)
    
    run_sep <- function(covars){
      mod <- get_form(spec_id, covars) %>%
        as.formula() %>%
        lm(df_m)
      
      coefs <- coef(mod)
      
      c(coefs["sep_ridit"],
        r2 = broom::glance(mod)[[1]]) %>%
        enframe(name = "term", value = "estimate")
    }
    
    r2_prs <- get_form(spec_id, "prs") %>%
      as.formula() %>%
      lm(df_m) %>%
      broom::glance(mod) %>%
      pull(1)
    
    bind_rows(bivar = run_sep("sep_ridit"),
              adjust = run_sep(c("prs", "sep_ridit")),
              .id = "mod") %>%
      uncount(ifelse(term == "r2", 2, 1), .id = "id") %>%
      mutate(term = ifelse(id == 2, "r2_diff", term),
             estimate = ifelse(id == 2, estimate - r2_prs, estimate))
  }
  
  
  map_dfr(1:500, get_boot, .id = "boot") %>%
    group_by(mod, term) %>%
    summarise(get_ci(estimate),
              .groups = "drop") %>%
    mutate(n = nrow(df_mod))
  
}

tic()
res_sep <- get_furrr(mod_specs, get_sep)
toc()


# 6. SEP Multiplicative
sep_vals <- count(df_long, sep_ridit) %>%
  drop_na() %>%
  pull(1) %>%
  c(0, ., 1)

get_mult <- function(spec_id){
  
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id)
  
  mod <- get_form(spec_id, "prs*sep_ridit") %>%
    as.formula() %>%
    lm(df_mod)
  
  mrg <- margins(mod, at = list(sep_ridit = sep_vals), variables = "prs") %>%
    tidy(conf.int = TRUE) %>%
    select(sep_ridit = 3, beta = 4, p = 7, lci = 8, uci = 9) %>%
    list()
  
  res <- tidy(mod, conf.int = TRUE) %>%
    filter(term == "prs:sep_ridit") %>%
    select(beta = 2, p = 5, lci = 6, uci = 7) %>%
    bind_cols(glance(mod) %>% select(r2 = 1, n = 12)) %>%
    list()
  
  tibble(margins = mrg, mod = res)
  
}

tic()
res_mult <- get_furrr(mod_specs, get_mult)
toc()

rm(sep_vals)

# 8. Save Objects ----
save(res_prs, res_attrit, res_quant,
     res_sep, res_mult, mod_specs, get_ci,
     file = "Data/regression_results.Rdata")
