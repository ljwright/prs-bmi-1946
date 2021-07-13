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
plan(multisession)

# 1. Load Data ----
load("Data/model_objects.Rdata")

# 2. Model Functions ----
get_spec <- function(spec_id){
  spec <- mod_specs %>%
    filter(spec_id == !!spec_id) %>%
    as.list()
  
  return(spec)
}

get_df <- function(spec_id){
  spec <- get_spec(spec_id)
  
  df_mod <- df_bmi %>%
    left_join(df_prs, by = "id") %>%
    filter(sex %in% sexes[[!!spec$sex]],
           age == !!spec$age) %>%
    rename(prs = all_of(spec$prs))
  
  if (spec$obs == "cc"){
    df_mod <- inner_join(df_mod, id_all, by = "id")
  }
  
  return(df_mod)
}

get_form <- function(spec_id, covars = "prs"){
  spec <- get_spec(spec_id)
  
  mod_form <- glue_collapse(covars, ' + ') %>%
    paste(spec$bmi_var, "~", .)
  
  if (spec$sex == "all") mod_form <- glue("{mod_form} + sex")
  
  return(mod_form)
}

get_ci <- function(estimates){
  quantile(estimates, c(.5, .025, .975)) %>%
    as_tibble_row() %>%
    rename(beta = 1, lci = 2, uci = 3)
}


# 3. Linear Regression Models ----
get_lm <- function(spec_id){
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id) %>%
    select(all_of(spec$bmi_var), prs, sex) %>%
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
df_lm <- mod_specs %>%
  mutate(res = future_map(spec_id, get_lm, 
                          .progress = TRUE,
                          .options = furrr_options(seed = NULL))) %>%
  unnest(res)
toc()


# 4. Quantile Regression ----
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

df_quant <- mod_specs %>%
  mutate(res = map(spec_id, get_quant)) %>%
  unnest(res)


# 5. SEP Additive ----
get_sep <- function(spec_id){
  
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id) %>% 
    select(all_of(spec$bmi_var), prs, fsc4_ridit, sex) %>%
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
      
      c(coefs["fsc4_ridit"],
        r2 = broom::glance(mod)[[1]]) %>%
        enframe(name = "term", value = "estimate")
    }
    
    r2_prs <- get_form(spec_id, "prs") %>%
      as.formula() %>%
      lm(df_m) %>%
      broom::glance(mod) %>%
      pull(1)
    
    bind_rows(bivar = run_sep("fsc4_ridit"),
              adjust = run_sep(c("prs", "fsc4_ridit")),
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
df_sep <- mod_specs %>%
  mutate(res = future_map(spec_id, get_sep, 
                          .progress = TRUE,
                          .options = furrr_options(seed = NULL))) %>%
  unnest(res)
toc()


# 6. SEP Multiplicative
fsc4_vals <- count(df_prs, fsc4_ridit) %>%
  drop_na() %>%
  pull(1) %>%
  c(0, ., 1)

get_mult <- function(spec_id){
  
  spec <- get_spec(spec_id)
  
  df_mod <- get_df(spec_id) %>% 
    select(all_of(spec$bmi_var), prs, fsc4_ridit, sex) %>%
    drop_na()
  
  mod <- get_form(spec_id, "prs*fsc4_ridit") %>%
    as.formula() %>%
    lm(df_mod)
  
  mrg <- margins(mod, at = list(fsc4_ridit = fsc4_vals), variables = "prs") %>%
    tidy(conf.int = TRUE) %>%
    select(fsc4_ridit = 3, beta = 4, p = 7, lci = 8, uci = 9) %>%
    list()
  
  res <- tidy(mod, conf.int = TRUE) %>%
    filter(term == "prs:fsc4_ridit") %>%
    select(beta = 2, p = 5, lci = 6, uci = 7) %>%
    bind_cols(glance(mod) %>% select(r2 = 1, n = 12)) %>%
    list()
  
  tibble(margins = mrg, mod = res)
  
}

df_mult <- mod_specs %>%
  mutate(res = map(spec_id, get_mult)) %>%
  unnest(res)

rm(fsc4_vals)


# 7. Splines ----
get_splines <- function(spec_id){
  df_mod <- get_df(spec_id)
  
  df_s <- df_mod %>%
    distinct(prs) %>%
    drop_na() %>%
    arrange(prs) %>%
    mutate(ns(prs, 2) %>%
             as_tibble() %>%
             rename_with(~ glue("ns_{.x}")))
  
  df_mod <- left_join(df_mod, df_s, by = "prs") %>%
    select(bmi, sex, matches("^ns_")) %>%
    drop_na()
  
  df_s <- df_s %>%
    filter(row_number() %% 10 == 1) %>%
    pivot_longer(-prs, names_to = "term")
  
  mod_form <- get_form(spec_id, str_subset(names(df_mod), "^ns_"))
  
  get_boot <- function(boot){
    set.seed(boot)
    df_m <- sample_frac(df_mod, replace = TRUE)
    
    as.formula(mod_form) %>%
      lm(df_m) %>%
      coef()
  }
  boots <- map(1:500, get_boot)
  
  list(boots = boots, splines = df_s)
}

tic()
df_splines <- mod_specs %>%
  filter(bmi_var == "bmi") %>%
  mutate(res = future_map(spec_id, get_splines, 
                          .progress = TRUE,
                          .options = furrr_options(seed = NULL)))
toc()

# Observed Data
get_splines_obs <- function(spec_id){
  df_mod <- get_df(spec_id) 
  mod_form <- get_form(spec_id)
  
  df_s <- df_mod %>%
    distinct(prs) %>%
    drop_na() %>%
    arrange(prs) %>%
    mutate(ns(prs, 2) %>%
             as_tibble() %>%
             rename_with(~ glue("ns_{.x}")))
  
  df_mod <- left_join(df_mod, df_s, by = "prs") %>%
    select(bmi, sex, matches("^ns_")) %>%
    drop_na()
  
  df_s <- df_s %>%
    filter(row_number() %% 10 == 1) %>%
    pivot_longer(-prs, names_to = "term")
  
  mod_form <- get_form(spec_id, str_subset(names(df_mod), "^ns_"))
  
  coefs <- as.formula(mod_form) %>%
    lm(df_mod) %>%
    coef()
  
  df_s %>%
    group_by(prs) %>%
    summarise(beta = sum(value*coefs[term]),
              .groups = "drop") %>%
    rename(prs_val = prs)
}

df_splines_ob <- mod_specs %>%
  filter(bmi_var == "bmi") %>%
  mutate(res = map(spec_id, get_splines_obs)) %>%
  unnest(res)


# 8. Save Objects ----
save(df_lm, df_quant,
     df_sep, df_mult,
     df_splines, df_splines_ob,
     mod_specs, get_ci,
     file = "Data/regression_results.Rdata")
