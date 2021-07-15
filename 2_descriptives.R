library(tidyverse)
library(haven)
library(glue)
library(broom)
library(gallimaufr)
library(magrittr)
library(corrr)
library(ggridges)

rm(list = ls())

# DECIDE KERNAL DENSITY AND SCATTER SAMPLE AND WAY OF MEASURING BMI
# DENSITY OF PRS BY FOLLOW-UP
# DENSITY OF BMI BY WHETHER PRS OBSERVED

# 1. Load Data ----
load("Data/df_long.Rdata")

df_long <- df_long %>%
  drop_na(matches("prs"), bmi) %>%
  count(id) %>%
  full_join(df_long, by = "id") %>%
  mutate(sample = if_else(n == length(unique(df_long$age)), "cc", "obs", "obs")) %>%
  select(-n)%>%
  mutate(age_f = ordered(age))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gwas_dict <- c(prs_k = "Khera et al. (2019)", 
               prs_r = "Richardson et al. (2020)", 
               prs_v = "Vogelezang et al. (2020)")

save(cbbPalette, gwas_dict, file = "Data/helpers.Rdata")


# 2. Kernel Density BMI ----
plot_bmi <- function(df){
  ggplot(df) +
    aes(x = bmi) +
    facet_wrap(~ age_f, ncol = 4) +
    geom_density(data = select(df_long, -age_f), aes(group = age),
                 color = "grey70", fill = "grey70", alpha = 0.4) +
    geom_density(aes(color = age_f, fill = age_f), alpha = 0.7) +
    guides(color = FALSE, fill = FALSE) +
    theme_bw() +
    labs(x = "BMI", y = "Density")
}

plot_bmi(df_long)
ggsave("Images/density_obs.png",
       height = 21, width = 29.7, units = "cm")

df_long %>%
  filter(sample == "cc") %>%
  plot_bmi() 
ggsave("Images/density_cc.png",
       height = 21, width = 29.7, units = "cm")

# 2. PRS, BMI, Height, Weight Correlations ----
# PRS-BMI Scatter Plot
df_scat <- df_long %>%
  select(age, matches("prs"), bmi) %>%
  pivot_longer(matches("prs"), names_to = "prs_var", values_to = "prs_value") %>%
  drop_na() %>%
  group_by(age, prs_var) %>%
  mutate(bmi = wtd_scale(bmi)) %>%
  ungroup()

plot_scatter <- function(prs_var){
  p <- df_scat %>%
    filter(prs_var == !!prs_var) %>%
    ggplot() +
    aes(x = prs_value, y = bmi) +
    facet_wrap(~ age) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_jitter(data = filter(df_scat, prs_var == "prs_k") %>% select(-age),
                alpha = 0.1, color = "grey80") +
    geom_jitter(alpha = 0.2, color = cbbPalette[4]) 
  
  glue("Images/scatter_bmi_{prs}.png") %>%
    ggsave(p, height = 21, width = 29.7, units = "cm")
  
  return(p)
}

rm(df_scat)

# PRS and BMI/Height/Weight Correlations
df_long %>%
  select(id, matches("prs_"),
         age, bmi, height, weight) %>%
  drop_na(bmi) %>%
  pivot_longer(c(bmi, height, weight), names_to = "phenotype", values_to = "pheno_value") %>%
  pivot_longer(matches("prs"), names_to = "prs", values_to = "prs_value") %>%
  drop_na() %>%
  mutate(pheno_clean = ifelse(phenotype == "bmi", 
                              str_to_upper(phenotype),
                              str_to_title(phenotype)),
         prs_clean = factor(gwas_dict[prs], gwas_dict)) %>%
  group_by(age, pheno_clean, prs_clean) %>%
  summarise(corr = cor(pheno_value, prs_value),
            .groups = "drop") %>%
  mutate(age = factor(age)) %>%
  ggplot() +
  aes(x = age, y = corr, color = pheno_clean,
      shape = pheno_clean, group = pheno_clean) +
  facet_grid(prs_clean ~ ., switch = "y") +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = cbbPalette[6:8]) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0)) +
  labs(x = "Age", y = "Correlation", 
       color = NULL, shape = NULL)
ggsave("Images/prs_corr.png",
       height = 16, width = 21, units = "cm")

# BMI and Height/Weight Correlations
df_long %>%
  select(age, bmi, weight, height) %>%
  drop_na() %>%
  arrange(age) %>%
  group_split(age, .keep = FALSE) %>%
  map_dfr(~ correlate(.x, quiet = TRUE) %>%
            stretch(),
          .id = "age") %>%
  drop_na() %>%
  mutate(age = unique(df_long$age)[as.numeric(age)] %>%
           as.factor(),
         phenotype = str_to_title(y)) %>%
  filter(x == "bmi") %>%
  ggplot() +
  aes(x = age, y = r, group = phenotype,
      color = phenotype, shape = phenotype) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = cbbPalette[7:8]) +
  scale_shape_manual(values = c(17, 15)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Age", y = "Correlation", 
       color = NULL, shape = NULL)
ggsave("Images/bmi_corr.png",
       height = 9.9, width = 21, units = "cm")

# 2. Attrition ----
# Distribution of BMI by PRS Missing/Observed
df_attrit <- df_long %>%
  select(id, age, matches("prs"), bmi) %>%
  pivot_longer(matches("prs"), names_to = "prs_var", values_to = "prs") %>%
  mutate(across(c(bmi, prs),
                list(miss = ~ if_else(is.na(.x), "Missing", "Observed") %>%
                  paste(str_to_upper(cur_column())) %>%
                  factor()),
                .names = "{.fn}_{.col}"))

plot_attrit <- function(x_var, prs = "prs_k"){
  miss_var <- ifelse(x_var == "bmi", "miss_prs", "miss_bmi")
  x_lab <- ifelse(x_var == "bmi", "BMI", "Polygenic Risk Score")

  p <- df_attrit %>%
    filter(prs_var == !!prs) %>%
    ggplot() +
    aes_string(x = x_var, color = miss_var, fill = miss_var) +
    facet_wrap(~ age, scales = "free") +
    geom_density(alpha = 0.3) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = x_lab, y = "Density",
         color = NULL, fill = NULL)
  
  glue("Images/attrit_{miss_var}_{prs}.png") %>%
    ggsave(p, height = 21, width = 29.7, units = "cm")
  
  return(p)
}

distinct(df_x, prs_var) %>%
  expand_grid(x_var = c("bmi", "prs")) %$%
  map2(x_var, prs_var, plot_attrit)

rm(df_attrit)


# 3. PRS x Socio-Economic Position -----
df_sep <- df_long %>%
  distinct(id, sep, across(matches("prs"))) %>%
  pivot_longer(matches("prs"), names_to = "prs", values_to = "prs_value") %>%
  mutate(gwas_clean = factor(gwas_dict[prs], gwas_dict),
         sep_f = sep) %>%
  drop_na()

df_sep %>%
  mutate(sep = fct_rev(sep)) %>%
  ggplot() +
  aes(x = prs_value, y = sep) +
  facet_grid(gwas_clean ~ ., switch = "y") +
  geom_density_ridges(fill = cbbPalette[6], alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0)) +
  labs(x = "Polygenic Risk Score", y = NULL)
ggsave("Images/prs_ridges_sep.png",
       height = 21, width = 21, units = "cm")

ggplot(df_sep) +
  aes(x = prs_value, group = sep) +
  facet_grid(gwas_clean ~ sep_f, switch = "y") +
  geom_density(data = select(df_sep, -sep_f),
               color = "grey70", fill = "grey70", alpha = 0.4) +
  geom_density(aes(color = sep, fill = sep), alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Polygenic Risk Score", y = NULL) +
  guides(color = FALSE, fill = FALSE)
ggsave("Images/prs_density_sep.png",
       height = 21, width = 29.7, units = "cm")

df_sep %>%
  nest(data = -gwas_clean) %>%
  mutate(res = map(data,
                   ~ lm(prs_value ~ sep, .x) %>%
                     tidy(conf.int = TRUE))) %>%
  unnest(res) %>%
  select(-data) %>%
  filter(str_detect(term, "^sep")) %>%
  mutate(term = str_replace(term, "^sep", "") %>%
           ordered(levels(df_sep$sep)) %>%
           fct_rev()) %>%
  ggplot() +
  aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high) +
  facet_wrap(~ gwas_clean) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_pointrange(color = cbbPalette[7]) +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, y = "Difference in Polygenic Risk Score")
ggsave("Images/mean_prs_sep.png",
       height = 16, width = 21, units = "cm")

df_sep %>%
  mutate(sep = fct_rev(sep)) %>%
  ggplot() +
  aes(x = sep, y = prs_value) +
  facet_wrap(~ gwas_clean) +
  geom_jitter(alpha = 0.3, color = cbbPalette[2]) +
  geom_boxplot(fill = NA, color = "grey50") +
  geom_violin(fill = NA, color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, y = "Polygenic Risk Score")
ggsave("Images/violin_prs_sep.png",
       height = 9.9, width = 21, units = "cm")
