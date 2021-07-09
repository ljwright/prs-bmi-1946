library(tidyverse)
library(glue)
library(magrittr)
library(paletteer)

rm(list = ls())

# 1. Load Data ----
load("Data/regression_results.Rdata")

gwas_dict <- c(zprs_k = "Khera et al. (2019)", 
               zprs_r = "Richardson et al. (2020)", 
               zprs_v = "Vogelezang et al. (2020)")

# 2. Linear Regressions ----
plot_lm <- function(obs, bmi_var, save_p = FALSE){
  p <- df_lm %>%
    filter(obs == !!obs, bmi_var == !!bmi_var) %>%
    mutate(age = factor(age) %>% fct_rev(),
           sex = str_to_title(sex),
           stat = ifelse(term == "prs", "Effect Size", "R-Squared"),
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

df_lm %>%
  distinct(obs, bmi_var) %$%
  map2(obs, bmi_var, plot_lm, TRUE)


# 3. Quantile Regressions ----
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

# 4. SEP Additive ----
mod_dict <- c(bivar = "Bivariate", adjust = "+ Polygenic Risk Score")
class_dict <- c("1" = "I", "2" = "II", "3" = "III", "4" = "IV", "v unskilled" = "V")

plot_sep <- function(obs, sex, bmi_var, save_p = FALSE){
  df_res <- df_sep %>%
    filter(sex == !!sex, obs == !!obs, bmi_var == !!bmi_var,
           !(term == "r2_diff" & mod == "bivar")) %>%
    mutate(age = factor(age),
           mod = ifelse(term == "r2_diff", NA, mod),
           mod_clean = factor(mod_dict[mod], mod_dict),
           sex = str_to_title(sex),
           gwas_clean = factor(gwas_dict[prs], gwas_dict))
  
  
  p1 <- df_res %>%
    filter(str_detect(term, "^fsc4")) %>%
    mutate(term = str_replace(term, "fsc4", ""),
           class_clean = factor(class_dict[term], class_dict, ordered = TRUE)) %>%
    ggplot() +
    aes(x = age, y = beta, ymin = lci, ymax = uci,
        color = class_clean, fill = class_clean, group = class_clean) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(gwas_clean ~ mod_clean, scales = "free_x", switch = "y") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom") +
    labs(x = "Age", y = NULL, color = "Social Class",
         fill = "Social Class")
  
  
  p2 <- df_res %>%
    filter(term == "r2") %>%
    mutate(term_clean = ifelse(term == "r2", "R-Squared", "Addtional Variance Explained") %>%
             fct_rev()) %>%
    ggplot() +
    aes(x = age, y = beta, ymin = lci, ymax = uci, color = mod_clean) +
    facet_grid(gwas_clean ~ term_clean, scales = "free_x") +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          strip.text.y  = element_blank()) +
    labs(x = "Age", y = NULL, color = "Model")
  
  p3 <- df_res %>%
    filter(term == "r2_diff") %>%
    mutate(term_clean = ifelse(term == "r2", "R-Squared", "Addtional Variance Explained") %>%
             fct_rev()) %>%
    ggplot() +
    aes(x = age, y = beta, ymin = lci, ymax = uci) +
    facet_grid(gwas_clean ~ term_clean, scales = "free_x") +
    geom_pointrange(color = paletteer_d("RColorBrewer::Set1")[3]) +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          strip.text.y  = element_blank()) +
    labs(x = "Age", y = NULL, color = "Model")
  
  p <- p1 + p2 + p3 +  plot_layout(widths = c(2, 1, 1))
  
  if (save_p){
    glue("Images/sep_{obs}_{bmi_var}_{sex}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

df_quant %>%
  distinct(obs, sex, bmi_var) %$%
  pmap(list(obs, sex, bmi_var), plot_sep, TRUE)


# 5. SEP Multiplicative ----
plot_quant <- function(obs, bmi_var, save_p = FALSE){
  p <- df_mult %>%
    filter(obs == !!obs, bmi_var == !!bmi_var) %>%
    mutate(age = factor(age),
           sex = str_to_title(sex),
           gwas_clean = factor(gwas_dict[prs], gwas_dict),
           term = str_replace(term, "prs\\:fsc4", ""),
           class_clean = factor(class_dict[term], class_dict, ordered = TRUE)) %>%
    ggplot() +
    aes(x = age, y = beta, ymin = lci, ymax = uci,
        color = class_clean, fill = class_clean, group = class_clean) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(gwas_clean ~ sex, scales = "free_x", switch = "y") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom") +
    labs(x = "Age", y = NULL, color = "Social Class",
         fill = "Social Class")
  
  if (save_p){
    glue("Images/mult_{obs}_{bmi_var}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

df_mult %>%
  distinct(obs, bmi_var) %$%
  map2(obs, bmi_var, plot_quant, TRUE)


# 6. Splines ----
make_splines <- function(res){
  map_dfr(res$boots,
          ~ enframe(.x, name = "term", value = "coef"),
          .id = "boot") %>%
    mutate(boot = as.integer(boot)) %>%
    right_join(res$splines, by = "term") %>%
    group_by(boot, prs) %>%
    summarise(estimate = sum(value*coef),
              .groups = "drop") %>%
    rename(prs_val = prs)
}

df_pred <- df_splines %>%
  mutate(res = map(res, make_splines))

df_ci <- df_pred %>%
  mutate(res = map(res, 
                   ~ .x %>%
                     group_by(prs_val) %>%
                     summarise(get_ci(estimate)))) %>%
  unnest(res)

plot_splines <- function(obs, prs, save_p = FALSE){
  p <- df_ci %>%
    filter(obs == !!obs, prs == !!prs) %>%
    mutate(age = factor(age) %>% ordered(),
           sex = str_to_title(sex),
           gwas_clean = factor(gwas_dict[prs], gwas_dict)) %>%
    ggplot() +
    aes(x = prs_val, y = beta, ymin = lci, ymax = uci,
        color = sex, fill = sex) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ age, scales = "free_y") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom") +
    labs(x = "Polygenic Risk Score", y = "Marginal Effect", 
         color = NULL, fill = NULL)
  
  if (save_p){
    glue("Images/splines_{obs}_{prs}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

df_ci %>%
  distinct(obs, prs) %$%
  map2(obs, prs, plot_splines, TRUE)

# Second Figure
plot_splines_2 <- function(obs, prs, save_p = FALSE){
  clean_res <- function(df_res, obs, prs){
    df_res %>%
      filter(obs == !!obs, prs == !!prs) %>%
      mutate(age = factor(age) %>% ordered(),
             sex = str_to_title(sex),
             gwas_clean = factor(gwas_dict[prs], gwas_dict))
  }
  
  p <- clean_res(df_pred, obs,prs) %>%
    unnest(res) %>%
    mutate(spec_id = spec_id + boot/1000) %>%
    ggplot() +
    aes(x = prs_val, y = estimate, color = sex, group = spec_id) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ age, scales = "free_y") +
    geom_line(alpha = 0.05) +
    geom_line(data = clean_res(rename(df_splines_ob, estimate = beta), obs, prs),
              size = 1.25) +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom") +
    labs(x = "Polygenic Risk Score", y = "Marginal Effect", 
         color = NULL, fill = NULL) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))
  
  if (save_p){
    glue("Images/splines2_{obs}_{prs}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

df_ci %>%
  distinct(obs, prs) %$%
  map2(obs, prs, plot_splines_2, TRUE)
