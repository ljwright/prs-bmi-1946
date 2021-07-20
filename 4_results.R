library(tidyverse)
library(glue)
library(magrittr)
library(paletteer)
library(patchwork)

rm(list = ls())

# 1. Load Data ----
load("Data/df_long.Rdata")
load("Data/regression_results.Rdata")
load("Data/helpers.Rdata")

dep_dict <- c(raw = "Raw BMI", std = "Standardised",
              ln = "% Change", rank = "Percentile Rank")

clean_res <- function(df_res){
  df_res %>%
    mutate(age_f = ordered(age),
           age_rev = fct_rev(age_f),
           sex_clean = str_to_title(sex),
           corrected = str_detect(dep_var, "_corrected"),
           dep_type = str_replace(dep_var, "_corrected", "") %>%
             str_replace("bmi_", ""),
           dep_clean = factor(dep_dict[dep_type], dep_dict),
           prs_clean = factor(gwas_dict[prs_var], gwas_dict))
}


# 2. Main Regressions ----
res_prs <- clean_res(res_prs) %>%
  mutate(term_clean = ifelse(term == "prs", "Effect Size", "R-Squared"))

plot_tbl <- df_long %>%
  select(age, matches("prs"), bmi, female) %>%
  pivot_longer(matches("prs"), names_to = "prs_var", values_to = "prs") %>%
  drop_na() %>%
  group_by(age, prs_var) %>%
  summarise(`Mean BMI` = mean(bmi), SD = sd(bmi),
            .groups = "drop") %>%
  pivot_longer(c(`Mean BMI`, `SD`), values_to = "beta", names_to = "stat") %>%
  mutate(age = factor(age) %>% fct_rev(),
         beta = round(beta, 2)) %>%
  ggplot() +
  aes(x = age, label = beta, y = stat) +
  facet_grid(prs_var ~ stat, scales = "free") +
  geom_text(size = 3) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        strip.text.y = element_blank(),
        panel.grid = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"))

# Plots
plot_main <- function(df, facet_form, colors, legend_pos){
  ggplot(df) +
    aes(x = age_rev, y = beta, ymin = lci, ymax = uci, 
        color = sex_clean, shape = sex_clean) +
    facet_grid(facet_form, scales = "free_x", switch = "y") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_color_manual(values = cbbPalette[colors]) +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = legend_pos,
          panel.spacing = unit(1, "lines")) +
    labs(x = NULL, y = NULL, color = NULL, shape = NULL)
}

# Main Results
p <- res_prs %>%
  filter(dep_var == "bmi_raw", sex == "all") %>%
  plot_main("prs_clean ~ term_clean", 1, "off")
p + plot_tbl + patchwork::plot_layout(widths = c(4, 1.1))
ggsave("Images/prs_main.png", height = 21, width = 21, units = "cm")
rm(p, plot_tbl)

# Main Results by Sex
res_prs %>%
  filter(dep_var == "bmi_raw", sex != "all") %>%
  plot_main("prs_clean ~ term_clean", 6:7, "bottom")
ggsave("Images/prs_sex.png", height = 21, width = 21, units = "cm")

# Main Results x Definition of BMI
res_prs %>%
  filter(corrected == FALSE, str_detect(dep_var, "bmi"),
         sex == "all", term == "prs") %>%
  plot_main("prs_clean ~ dep_clean", 1, "off")
ggsave("Images/prs_all_dep.png", height = 21, width = 29.7, units = "cm")

# Main Results x Definition of Corrected BMI
res_prs %>%
  filter(corrected == TRUE, str_detect(dep_var, "bmi"),
         sex == "all", term == "prs") %>%
  plot_main("prs_clean ~ dep_clean", 1, "off")
ggsave("Images/prs_all_corrected.png", height = 21, width = 29.7, units = "cm")

# Main Results - Fat Mass, etc.
alt_dict <- c(fat_ratio = "Fat Ratio", fat_mass = "Fat Mass", lean_mass = "Lean Mass")

res_prs %>%
  filter(dep_var %in% c("fat_ratio", "fat_mass", "lean_mass"),
         sex == "all", term == "prs") %>%
  mutate(dep_clean = factor(alt_dict[dep_var], alt_dict)) %>%
  plot_main("~ prs_clean", 1, "off") +
  aes(x = dep_clean) +
  labs(x = NULL)
ggsave("Images/prs_alternative.png", height = 21, width = 29.7, units = "cm")

# 3. Attrition Regressions ----
res_attrit <- clean_res(res_attrit) %>%
  mutate(age_from = factor(age_from) %>%
           fct_recode("Observed" = "0", "Complete Cases" = "2"))

plot_attrit <- function(prs_var, dep_var, sex, save_p = FALSE){
  df_attrit <- res_attrit %>%
    filter(dep_var == !!dep_var, 
           prs_var == !!prs_var,
           sex == !!sex, 
           age_from != "69")
  
  p <- ggplot(df_attrit) +
    aes(x = age_f, y = beta, ymin = lci, ymax = uci) +
    facet_wrap(~ age_from) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(aes(group = age_from), alpha = 0.2,
                fill = cbbPalette[6]) +
    geom_line(data = rename(df_attrit, age_samp = age_from),
              aes(group = age_samp), color = "grey70") +
    geom_line(aes(group = age_from), 
              color = cbbPalette[6])  +
    theme_bw() +
    theme(panel.spacing = unit(1, "lines")) +
    labs(x = "Age", y = "Marginal Effect")
  
  if (save_p){
    glue("Images/attrit_{prs_var}_{sex}_{dep_var}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

res_attrit %>%
  distinct(prs_var, dep_var, sex) %>%
  filter(sex == "all",
         dep_var %in% c("bmi_raw", "bmi_std")) %$%
  pwalk(list(prs_var, dep_var, sex), plot_attrit, TRUE)


# 4. Quantile Regressions ----
res_quant <- clean_res(res_quant) %>%
  mutate(tau_clean = glue("{tau*100}th") %>%
           fct_reorder(tau))

plot_quant <- function(dep_var, sex, prs_var, save_p = FALSE){
  df_res <- res_quant %>%
    filter(sex == !!sex, dep_var == !!dep_var, prs_var == !!prs_var)
  
  p <- ggplot(df_res) +
    aes(x = tau_clean, y = beta, ymin = lci, ymax = uci, group = age) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ age_f) +
    geom_line(data = select(df_res, -age_f), color = "grey70") +
    geom_ribbon(color = NA, fill = cbbPalette[6], alpha = 0.2) +
    geom_line(color = cbbPalette[6]) +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing = unit(1, "lines")) +
    labs(x = "Percentile", y = "Marginal Effect",
         color = "Age", fill = "Age") +
    guides(color = guide_legend(nrow = 2), 
           fill = guide_legend(nrow = 2))
  
  if (save_p){
    glue("Images/quant_{sex}_{dep_var}_{prs_var}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

res_quant %>%
  distinct(dep_var, sex, prs_var) %>%
  filter(sex == "all", dep_var == "bmi_raw") %$%
  pwalk(list(dep_var, sex, prs_var), plot_quant, TRUE)


# 5. SEP Additive ----
mod_dict <- c(bivar = "Bivariate", adjust = "+ Polygenic Risk Score")

sep_dict <- c(sep_ridit = "Father's Social Class (Ridit)",
              medu = "Mother's Education",
              medu_ridit = "Mother's Education (Ridit)")

res_sep <- clean_res(res_sep) %>%
  mutate(mod_clean = factor(mod_dict[mod], mod_dict),
         sep_clean = factor(sep_dict[sep_var], sep_dict))

res_sep %>%
  filter(sex == "all", dep_var == "bmi_raw") %>%
  distinct(age_f, prs_clean, mod_clean, r_diff, sep_clean) %>%
  ggplot() +
  aes(x = age_f, y = r_diff, color = mod_clean, 
      group = mod_clean, shape = mod_clean) +
  facet_grid(prs_clean ~ sep_clean, scales = "free", switch = "y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  scale_color_manual(values = cbbPalette[6:7]) +
  scale_fill_manual(values = cbbPalette[6:7]) +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines")) +
  labs(x = "Age", y = NULL, color = "Model", shape = "Model")
ggsave("Images/sep_r.png", height = 21, width = 29.7, units = "cm")

plot_sep <- function(dep_var, sex, save_p = FALSE){
  p <- res_sep %>%
    filter(dep_var == !!dep_var, sex == !!sex) %>%
    ggplot() +
    aes(x = age_f, y = beta, ymin = lci, ymax = uci, group = mod_clean, 
        color = mod_clean, fill = mod_clean, linetype = mod_clean) +
    facet_grid(prs_clean ~ sep_clean, scales = "free", switch = "y") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    scale_color_manual(values = cbbPalette[6:7]) +
    scale_fill_manual(values = cbbPalette[6:7]) +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          panel.spacing = unit(1, "lines")) +
    labs(x = "Age", y = NULL, color = "Model", 
         fill = "Model", linetype = "Model")
  
  if (save_p){
    glue("Images/sep_{sex}_{dep_var}.png") %>%
      ggsave(p, height = 16, width = 21, units = "cm")
  }
  
  return(p)
}

res_sep %>%
  distinct(dep_var, sex) %>%
  filter(sex == "all",
         dep_var %in% c("bmi_raw", "bmi_std")) %$%
  walk2(dep_var, sex, plot_sep, TRUE)


# 6. SEP Multiplicative ----
res_mult <- res_mult %>%
  pivot_longer(c(margins, mod)) %>%
  unnest(value) %>%
  mutate(sep_level = case_when(sep_values == 0 ~ "Lowest SEP",
                               sep_values == 1 ~ "Highest SEP",
                               is.na(sep_values) ~ "Interaction Term",
                               TRUE ~ NA_character_) %>%
           factor(c("Lowest SEP", "Highest SEP", "Interaction Term")),
         sep_clean = factor(sep_dict[sep_var], sep_dict)) %>%
  filter(!is.na(sep_level)) %>%
  clean_res()

plot_mult <- function(dep_var, sex, save_p = FALSE){
  p <- res_mult %>%
    filter(dep_var == !!dep_var, sex == !!sex,
           is.na(sep_values)) %>%
    ggplot() +
    aes(x = age_f, y = beta, ymin = lci, ymax = uci, group = sep_clean) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(prs_clean ~ sep_clean, scales = "free_x", switch = "y") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          panel.spacing = unit(1, "lines")) +
    labs(x = "Age", y = NULL, color = NULL,
         fill = NULL)
  
  if (save_p){
    glue("Images/mult_{sex}_{dep_var}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

res_mult %>%
  distinct(dep_var, sex) %>%
  filter(sex == "all",
         dep_var %in% c("bmi_raw", "bmi_std")) %$%
  walk2(dep_var, sex, plot_mult, TRUE)
