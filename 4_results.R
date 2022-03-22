library(tidyverse)
library(glue)
library(magrittr)
library(paletteer)
library(patchwork)
library(gallimaufr)
library(flextable)
library(officer)

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

# Results Table
all_dict <- c(bmi_raw = "Absolute BMI", bmi_std = "Standardized BMI",
              bmi_ln = "Log BMI", bmi_rank = "BMI Percentile",
              bmi_corrected_raw = "Absolute BMI (Corrected)", 
              bmi_corrected_std = "Standardized BMI (Corrected)",
              bmi_corrected_ln = "Log BMI (Corrected)", 
              bmi_corrected_rank = "BMI Percentile (Corrected)",
              fat_ratio = "Fat Ratio", fat_mass = "Fat Mass",
              lean_mass = "Lean Mass", weight = "Weight",
              height = "Height")

header_lbls <- list(dep_clean = "Outcome", age = "Age", 
                    all_prs_k = "Khera et al. (2019)",
                    all_prs_ksig = "Khera et al. (2019)*",
                    all_prs_r = "Richardson et al. (2020)",
                    all_prs_v = "Vogelezang et al. (2020)", 
                    female_prs_k = "Khera et al. (2019)",
                    female_prs_ksig = "Khera et al. (2019)*",
                    female_prs_r = "Richardson et al. (2020)",
                    female_prs_v = "Vogelezang et al. (2020)", 
                    male_prs_k = "Khera et al. (2019)",
                    male_prs_ksig = "Khera et al. (2019)*",
                    male_prs_r = "Richardson et al. (2020)",
                    male_prs_v = "Vogelezang et al. (2020)")

span_lbls <- list(dep_clean = "", age = "", 
                  all_prs_k = "All", all_prs_ksig = "All",
                  all_prs_r = "All", all_prs_v = "All", 
                  female_prs_k = "Female", female_prs_ksig = "Female",
                  female_prs_r = "Female", female_prs_v = "Female", 
                  male_prs_k = "Male", male_prs_ksig = "Male",
                  male_prs_r = "Male", male_prs_v = "Male")

get_flx <- function(term){
  res_prs %>%
    mutate(dep_clean = factor(all_dict[dep_var], all_dict),
           across(beta:uci, round, 2),
           string = glue("{beta} ({lci}, {uci})")) %>%
    filter(term == !!term) %>%
    arrange(dep_clean, age, sex, prs_var) %>%
    select(dep_clean, age, sex, prs_var, string) %>%
    pivot_wider(names_from = sex, values_from = string) %>%
    pivot_wider(names_from = prs_var, values_from = all:male) %>%
    make_flx(header_lbls, span_lbls) %>%
    vline(j = c(2, 6, 10), part = "body",
          border = fp_border(color="grey50", width = 1, style = "dashed"))
}

flx_prs <- get_flx("prs")
flx_prs
save_as_docx(flx_prs, path = "Tables/results_table.docx")

flx_r2 <- get_flx("r2")
flx_r2
save_as_docx(flx_r2, path = "Tables/r2_table.docx")

rm(flx_prs, flx_r2)


# Small Tables
make_tbl <- function(df_res, pivot_var){
  df_res %>%
    mutate(across(beta:uci, round, 2),
           string = glue("{beta} ({lci}, {uci})")) %>%
    select(age, var = all_of(!!pivot_var), string) %>%
    arrange(age, var) %>%
    pivot_wider(names_from = var, values_from = string)
}

res_prs %>%
  filter(sex == "all", term == "prs", str_detect(dep_var, "(fat|lean)")) %>%
  mutate(across(beta:uci, round, 2),
         string = glue("b = {beta}; 95% CI = {lci}, {uci}")) %>%
  select(dep_var, prs_var, string) %>%
  arrange(dep_var, prs_var) %>%
  pivot_wider(names_from = prs_var, values_from = string)


# Mean and SD Plot
plot_tbl <- function(prs_vars){
  df_long %>%
    select(age, matches("prs"), bmi, female) %>%
    pivot_longer(matches("prs"), names_to = "prs_var", values_to = "prs") %>%
    drop_na() %>%
    filter(prs_var %in% prs_vars) %>%
    group_by(age, prs_var) %>%
    summarise(`Mean BMI` = mean(bmi), SD = sd(bmi),
              .groups = "drop") %>%
    pivot_longer(c(`Mean BMI`, `SD`), values_to = "beta", names_to = "stat") %>%
    mutate(age = factor(age) %>% fct_rev(),
           beta = round(beta, 1)) %>%
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
}


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

# Main Results (Khera and Vogelezang)
p <- res_prs %>%
  filter(dep_var == "bmi_raw", sex == "all",
         prs_var %in% c("prs_k", "prs_v")) %>%
  plot_main("prs_clean ~ term_clean", 1, "off")
p + plot_tbl(c("prs_k", "prs_v")) + patchwork::plot_layout(widths = c(4, 1.1))
ggsave("Images/prs_main.png", height = 16, width = 21, units = "cm")

# Full
p <- res_prs %>%
  filter(dep_var == "bmi_raw", sex == "all") %>%
  plot_main("prs_clean ~ term_clean", 1, "off")
p + plot_tbl(names(gwas_dict)) + patchwork::plot_layout(widths = c(4, 1.1))
ggsave("Images/prs_all.png", height = 29.7, width = 21, units = "cm")
rm(p, plot_tbl)

# Main Results by Sex
res_prs %>%
  filter(dep_var == "bmi_raw", sex != "all") %>%
  plot_main("prs_clean ~ term_clean", 6:7, "bottom")
ggsave("Images/prs_sex.png", height = 29.7, width = 21, units = "cm")

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
  facet_grid(~ prs_clean) +
  aes(x = dep_clean) +
  labs(x = NULL)
ggsave("Images/prs_alternative.png", 
       height = 9.9, width = 16, units = "cm")

# 3. Attrition Regressions ----
res_attrit <- clean_res(res_attrit) %>%
  mutate(age_from = factor(age_from) %>%
           fct_recode("Observed" = "0", "Complete Cases" = "2")) %>%
  filter(age_from != "6")

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


# Attrition R2
plot_attrit_r2 <- function(prs_var, save_p = FALSE){
  df_attrit <- clean_res(res_attrit_r2) %>%
    filter(!(age_from %in% c(6, 69))) %>%
    mutate(age_from = factor(age_from) %>%
             fct_recode("Observed" = "0", "Complete Cases" = "2")) %>%
    filter(prs_var == !!prs_var)
  
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
    labs(x = "Age", y = bquote("Incremental R"^2))
  
  if (save_p){
    glue("Images/attrit_r2_{prs_var}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
  }
  
  return(p)
}

unique(res_attrit_r2$prs_var) %>%
  walk(plot_attrit_r2, TRUE)


# 4. Quantile Regressions ----
res_quant <- clean_res(res_quant) %>%
  mutate(tau_clean = glue("{tau*100}th") %>%
           fct_reorder(tau))

flx_q <- res_quant %>%
  filter(sex == "all", dep_var == "bmi_raw", prs_var == "prs_k") %>%
  make_tbl("tau_clean") %>%
  make_flx(list(age = "Age"))
flx_q
save_as_docx(flx_q, path = "Tables/quantile_results.docx")

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

# Heat Maps
plot_heat <- function(dep_var, sex, prs_var, save_p = FALSE){
  df_res <- res_quant %>%
    filter(sex == !!sex, dep_var == !!dep_var, prs_var == !!prs_var)
  
  p <- ggplot(df_res) +
    aes(x = tau_clean, y = age_f, color = beta, fill = beta) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0,2)) +
    scale_color_viridis_c(limits = c(0, 2)) +
    theme_minimal() +
    labs(x = "Quantile", y = "Age", 
         color = "Estimate", fill = "Estimate") +
    theme(legend.position = "top") +
    guides(fill = guide_colorbar(title.position = 'top', 
                                 title.hjust = .5,                                
                                 barwidth = unit(20, 'lines'), 
                                 barheight = unit(.5, 'lines')))
  
  if (save_p){
    glue("Images/heat_{sex}_{dep_var}_{prs_var}.png") %>%
      ggsave(p, height = 9.9, width = 16, units = "cm")
  }
  
  return(p)
}
res_quant %>%
  distinct(dep_var, sex, prs_var) %>%
  filter(sex == "all", dep_var == "bmi_raw") %$%
  pwalk(list(dep_var, sex, prs_var), plot_heat, TRUE)


# 5. SEP Additive ----
mod_dict <- c(bivar = "Bivariate", adjust = "+ Polygenic Risk Score")

sep_dict <- c(
  all = "Multiple Adjusted",
  father_class = "Father's Social Class",
  father_class_ridit = "Father's Social Class (Ridit)",
  father_edu_level = "Father's Education Level",
  father_edu_level_ridit = "Father's Education Level (Ridit)",
  father_edu_years = "Father's Education Years",
  father_edu_years_ridit = "Father's Education Years (Ridit)",
  mother_edu_level = "Mother's Education Level",
  mother_edu_level_ridit = "Mother's Education Level (Ridit)",
  mother_edu_years = "Mother's Education Years",
  mother_edu_years_ridit = "Mother's Education Years (Ridit)",
  own_class = "Social Class @ Age 53",
  own_class_ridit = "Social Class @ Age 53 (Ridit)",
  own_edu_level = "Own Education Level",
  own_edu_level_ridit = "Own Education Level (Ridit)"
)

res_sep <- clean_res(res_sep) %>%
  mutate(sep_var = map_chr(sep_var,
                           ~ ifelse(length(.x) == 1, .x, "all")),
         mod_clean = factor(mod_dict[mod], mod_dict),
         sep_clean = factor(sep_dict[sep_var], sep_dict))

# Table
sep_tbl <- res_sep %>%
  filter(sex == "all", dep_var == "bmi_raw", mod == "adjust") %>%
  distinct(prs_var, age, sep_var, r_prs_diff) %>%
  mutate(r_diff = glue("{round(100*r_prs_diff, 2)}%")) %>%
  select(prs_var, age, sep_var, r_diff) %>%
  arrange(prs_var, age, age, sep_var) %>%
  pivot_wider(names_from = sep_var, values_from = r_diff) %>%
  mutate(prs_var = gwas_dict[prs_var])

header_lbls <- tibble(var = names(sep_tbl),
                      name = case_when(var == "prs_var" ~ "Polygenic Score",
                                       var == "age" ~ "Age",
                                       TRUE ~ sep_dict[var])) %>%
  deframe() %>%
  as.list()

flx_sep <- make_flx(sep_tbl, header_lbls) %>%
  align(align = "center", part = "all")
flx_sep
save_as_docx(flx_sep, path = "Tables/sep_results.docx")

# R2 Plot
res_sep %>%
  filter(sex == "all", dep_var == "bmi_raw") %>%
  filter(sep_var %in% c("father_class_ridit", "mother_edu_years", "mother_edu_years_ridit")) %>%
  distinct(age_f, prs_clean, mod_clean, r_prs_diff, sep_clean) %>%
  ggplot() +
  aes(x = age_f, y = r_prs_diff, color = mod_clean, 
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

res_sep %>%
  filter(sex == "all", dep_var == "bmi_raw", 
         prs_var == "prs_k", mod == "adjust") %>%
  filter(!str_detect(sep_var, "(edu_years|class_ridit|level_ridit)$")) %>%
  distinct(age_f, prs_clean, mod_clean, r_prs_diff, sep_clean) %>%
  ggplot() +
  aes(x = age_f, y = r_prs_diff, group = mod_clean) +
  facet_wrap(~ sep_clean, ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.spacing = unit(1, "lines")) +
  labs(x = "Age", y = NULL)
ggsave("Images/sep_additive_r2.png", 
       height = 21, width = 29.7, units = "cm")

# Main Plot
plot_sep <- function(dep_var, sex, save_p = FALSE){
  p <- res_sep %>%
    filter(dep_var == !!dep_var, sex == !!sex, 
           prs_var == "prs_k",
           str_detect(sep_var, "ridit")) %>%
    ggplot() +
    aes(x = age_f, y = beta, ymin = lci, ymax = uci, group = mod_clean, 
        color = mod_clean, fill = mod_clean, linetype = mod_clean) +
    facet_wrap(~ sep_clean, labeller = label_wrap_gen(25), nrow = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    scale_color_manual(values = cbbPalette[6:7]) +
    scale_fill_manual(values = cbbPalette[6:7]) +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = "bottom",
          panel.spacing = unit(1, "lines")) +
    labs(x = "Age", y = NULL, color = "Model", 
         fill = "Model", linetype = "Model")
  
  if (save_p){
    glue("Images/sep_{sex}_{dep_var}.png") %>%
      ggsave(p, height = 21, width = 29.7, units = "cm")
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
    filter(!str_detect(sep_var, "edu_years")) %>%
    ggplot() +
    aes(x = age_f, y = beta, ymin = lci, ymax = uci, group = sep_clean) +
    facet_grid(sep_clean ~ prs_clean, scales = "free_x", switch = "y",
               labeller = label_wrap_gen(16)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
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


# Additive and Multiplicative
plot_both <- function(sep_vars, legend.position){
  bind_rows("Main Effects" = res_sep, 
            "Interaction Term" = filter(res_mult, name == "mod"),
            .id = "facet") %>%
    filter(sex == "all", prs_var == "prs_k",
           dep_var == "bmi_raw", sep_var %in% !!sep_vars) %>%
    mutate(facet = fct_rev(facet)) %>%
    ggplot() +
    aes(x = age_f, y = beta, ymin = lci, ymax = uci,
        group = mod_clean, color = mod_clean, fill = mod_clean) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(facet ~ ., scales = "free", switch = "y") +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    scale_color_manual(values = cbbPalette[6:7],
                       breaks = levels(res_sep$mod_clean)) +
    scale_fill_manual(values = cbbPalette[6:7],
                      breaks = levels(res_sep$mod_clean)) +
    theme_minimal() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          legend.position = legend.position,
          panel.spacing = unit(1, "lines")) +
    labs(x = "Age", y = NULL, 
         color = "Model", fill = "Model")
}

plot_both("father_class_ridit", c(0.15, 0.85))
ggsave("Images/both_father_class_ridit.png", 
       height = 16, width = 21, units = "cm")

names(sep_dict) %>%
  str_subset("ridit") %>%
  str_subset("father_class_ridit", TRUE) %>%
  plot_both("bottom") +
  facet_grid(facet ~ sep_clean, scales = "free", switch = "y",
             labeller = label_wrap_gen(20))
ggsave("Images/both_all_ridit.png", 
       height = 21, width = 29.7, units = "cm")

