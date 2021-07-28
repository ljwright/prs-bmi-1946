library(tidyverse)
library(haven)
library(glue)
library(broom.mixed)
library(lme4)
library(scales)
library(splines)

rm(list = ls())

load("Data/df_long.Rdata")

get_splines <- function(x, deg_free, stub = "ns"){
  ns(x, deg_free) %>%
    as_tibble() %>%
    mutate(across(everything(), as.numeric)) %>%
    rename_with(~ glue("{stub}_{.x}"))
}

df_reg <- df_long %>%
  drop_na(bmi, matches("prs")) %>%
  # select(id, age, bmi, matches("prs"), sep_ridit) %>%
  arrange(id, age) %>%
  mutate(age_s = rescale(age, from = c(0, 69)),
         get_splines(age_s, 2),
         ns_0 = 1)

mod <- lmer(bmi ~ -1 + ns_0 + ns_1 + ns_2 +
              (-1 + ns_0 + ns_1 + ns_2 | id),
            df_long, REML = TRUE)


ns_long <- df_long %>%
  select(age, matches("^ns_")) %>%
  distinct() %>%
  pivot_longer(-age, names_to = "term")


# Predicted Values 
df_re <- ranef(mod) %>%
  pluck("id") %>%
  as_tibble() %>%
  bind_cols(id = mod@frame$id %>% unique()) %>%
  pivot_longer(-id, names_to = "term", values_to = "re")

df_fe <- fixef(mod) %>%
  enframe(name = "term", value = "fe")


df_pred <- left_join(df_re, df_fe, by = "term") %>%
  left_join(ns_long, by = "term") %>%
  mutate(effect = value*(re+fe)) %>%
  group_by(age, id) %>%
  summarise(pred = sum(effect), .groups = "drop")

df_fixed <- left_join(df_fe, ns_long, by = "term") %>%
  group_by(age) %>%
  summarise(pred = sum(fe*value), .groups = "drop") %>%
  mutate(id = -1)

df_prs <- df_long %>%
  select(id, prs_k) %>%
  mutate(prs = cut_number(prs_k, 5) %>%
           as.numeric() %>% ordered()) %>%
  left_join(df_pred, by = "id")


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(df_prs) +
  aes(x = age, y = pred, group = id) +
  facet_wrap(~ prs) +
  geom_line(data = select(df_prs, -prs), alpha = 0.05, color = "grey70") +
  geom_line(aes(color = prs), alpha = 0.2) + # color = "grey70") +
  geom_line(data = df_fixed, color = "red") + #cbbPalette[6])
  theme_bw()

get_ci <- function(estimates){
  quantile(estimates, c(.5, .025, .975)) %>%
    as_tibble_row() %>%
    rename(beta = 1, lci = 2, uci = 3)
}

df_prs %>%
  group_by(prs, age) %>%
  summarise(get_ci(pred), .groups = "drop") %>%
  ggplot() +
  aes(x = age, y = beta, ymin = lci, ymax = uci,
      color = prs, fill = prs) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line()


# 3. Interaction Terms ----
mod_form <- "bmi ~ -1 + ns_0 + ns_1 + ns_2 + sep_var*prs_var +
  (-1 + ns_0 + ns_1 + ns_2 | id)"

get_int <- function(sep_var, prs_var){
  df_mod <- df_reg %>%
    rename(prs_var = all_of(!!prs_var),
           sep_var = all_of(!!sep_var))
  
  mod_form %>%
    as.formula() %>%
    lmer(df_mod) %>%
    tidy(conf.int = TRUE) %>%
    filter(term == "sep_var:prs_var") %>%
    select(beta = estimate, lci = conf.low, uci = conf.high)
}

df_reg %>%
  count(sep, sep_ridit)

expand_grid(sep_var = c("sep_ridit", "medu", "medu_ridit"),
            prs_var = c("prs_k", "prs_v", "prs_r")) %>%
  mutate(map2_dfr(sep_var, prs_var, get_int)) %>%
  ggplot() +
  aes(x = prs_var, y = beta, ymin = lci, ymax = uci) +
  facet_wrap(~ sep_var) +
  geom_hline(yintercept = 0) +
  geom_pointrange()


df_long %>%
  select(fat_mass, lean_mass) %>%
  drop_na() %>%
  summarise(cor = cor(fat_mass, lean_mass))
