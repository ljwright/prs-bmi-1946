library(tidyverse)
library(glue)
library(magrittr)
library(paletteer)
library(patchwork)

rm(list = ls())

# 1. Load Data ----
load("Data/regression_results.Rdata")
load("Data/helpers.Rdata")

# 2. Main Regressions ----
# CHANGE TO LINE AND RIBBON!
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
          legend.position = "bottom",
          panel.spacing = unit(1, "lines")) +
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


# 3. Attrition Regressions ----
df_res <- res_attrit %>%
  filter(prs_var == "prs_k", dep_var == "bmi", sex == "all", age_from < 69) %>%
  mutate(age_from = ordered(age_from))

ggplot(df_res) +
  aes(x = age, y = beta, ymin = lci, ymax = uci) +
  facet_wrap(~ age_from) +
  geom_ribbon(aes(fill = age_from), color = NA, alpha = 0.2) +
  geom_line(data = rename(df_res, age_f = age_from),
            aes(group = age_f)) +
  geom_line(aes(color = age_from)) 


# 4. Quantile Regressions ----
plot_quant <- function(obs, bmi_var, save_p = FALSE){
  p <- df_quant %>%
    filter(obs == !!obs, bmi_var == !!bmi_var) %>%
    mutate(age = factor(age) %>% ordered(),
           sex = str_to_title(sex),
           gwas_clean = factor(gwas_dict[prs], gwas_dict),
           tau = glue("{tau*100}th") %>%
             fct_reorder(tau)) %>%
    ggplot() +
    aes(x = tau, y = beta, ymin = lci, ymax = uci,
        color = age, fill = age, group = age) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(gwas_clean ~ sex, scales = "free", switch = "y") +
    geom_ribbon(color = NA, alpha = 0.1) +
    geom_line() +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing = unit(1, "lines")) +
    labs(x = "Percentile", y = NULL, color = "Age", fill = "Age") +
    guides(color = guide_legend(nrow = 2),
           fill = guide_legend(nrow = 2))
  
  if (save_p){
    glue("Images/quant_{obs}_{bmi_var}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

df_quant %>%
  distinct(obs, bmi_var) %$%
  map2(obs, bmi_var, plot_quant, TRUE)

# 5. SEP Additive ----
mod_dict <- c(bivar = "Bivariate", adjust = "+ Polygenic Risk Score")

plot_sep <- function(obs, sex, bmi_var, save_p = FALSE){
  df_res <- df_sep %>%
    filter(sex == !!sex, obs == !!obs, bmi_var == !!bmi_var,
           !(term == "r2_diff" & mod == "bivar")) %>%
    mutate(age = factor(age),
           mod = ifelse(term == "r2_diff", NA, mod),
           mod_clean = factor(mod_dict[mod], mod_dict),
           sex = str_to_title(sex),
           gwas_clean = factor(gwas_dict[prs], gwas_dict),
           term_clean = case_when(term == "r2" ~ "R-Squared",
                                  term == "r2_diff" ~ "Addtional Variance Explained by SEP",
                                  TRUE ~ "Marginal Effect") %>%
             fct_rev())
  
  p1 <- df_res %>%
    filter(str_detect(term, "^fsc4")) %>%
    ggplot() +
    aes(x = age, y = beta, ymin = lci, ymax = uci, 
        color = mod_clean, fill = mod_clean, group = mod_clean) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(gwas_clean ~ term_clean, scales = "free", switch = "y") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    scale_color_manual(values = cbbPalette[6:7]) +
    scale_fill_manual(values = cbbPalette[6:7]) +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          panel.spacing = unit(1, "lines")) +
    labs(x = "Age", y = NULL, color = "Model", fill = "Model")
  
  p2 <- df_res %>%
    filter(term == "r2") %>%
    ggplot() +
    aes(x = age, y = beta, ymin = lci, ymax = uci) +
    facet_grid(gwas_clean ~ term_clean) +
    geom_pointrange(aes(color = mod_clean), position = position_dodge(0.5)) +
    geom_pointrange(data = filter(df_res, term == "r2_diff"),
                    color = cbbPalette[8]) +
    scale_color_manual(values = cbbPalette[6:7]) +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          strip.text.y  = element_blank(),
          panel.spacing = unit(1, "lines")) +
    guides(color = FALSE) +
    labs(x = "Age", y = NULL) 
  
  p <- p1 + p2  + plot_layout(widths = c(1, 2)) +
    plot_layout(guides = "collect") & 
    theme(legend.position = 'bottom')
  
  if (save_p){
    glue("Images/sep_{obs}_{bmi_var}_{sex}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

df_sep %>%
  distinct(obs, sex, bmi_var) %$%
  pmap(list(obs, sex, bmi_var), plot_sep, TRUE)


# 6. SEP Multiplicative ----
plot_mult <- function(obs, bmi_var, save_p = FALSE){
  p <- df_mult %>%
    pivot_longer(c(margins, mod)) %>%
    unnest(value) %>%
    mutate(ridit = case_when(fsc4_ridit == 0 ~ "Lowest SEP",
                             fsc4_ridit == 1 ~ "Highest SEP",
                             is.na(fsc4_ridit) ~ "Interaction Term",
                             TRUE ~ NA_character_) %>%
             factor(c("Lowest SEP", "Highest SEP", "Interaction Term")),
           age = factor(age)) %>%
    filter(!is.na(ridit)) %>%
    filter(obs == !!obs, bmi_var == !!bmi_var) %>%
    mutate(sex = str_to_title(sex),
           # age = factor(age),
           gwas_clean = factor(gwas_dict[prs], gwas_dict)) %>%
    ggplot() +
    aes(x = age, y = beta, ymin = lci, ymax = uci,
        color = ridit, fill = ridit, group = ridit) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(gwas_clean ~ sex, scales = "free_x", switch = "y") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          panel.spacing = unit(1, "lines")) +
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
  map2(obs, bmi_var, plot_mult, TRUE)
