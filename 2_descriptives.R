library(tidyverse)
library(haven)
library(glue)
library(broom)
library(gallimaufr)
library(magrittr)
library(corrr)
library(ggridges)
library(summarytools)
library(officer)
library(flextable)

rm(list = ls())

# 1. Load Data ----
load("Data/df_long.Rdata")

df_desc <- df_long %>%
  drop_na(matches("prs"), bmi) %>%
  count(id) %>%
  full_join(df_long, by = "id") %>%
  mutate(sample = if_else(n == length(unique(df_long$age)), "cc", "obs", "obs")) %>%
  select(-n)%>%
  mutate(age_f = ordered(age)) %>%
  arrange(id, age)

rm(df_long)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gwas_dict <- c(prs_k = "Khera et al. (2019) Adult",  
               prs_r_adult = "Richardson et al. (2020) Adult",
               prs_v = "Vogelezang et al. (2020) Child",
               prs_r_child = "Richardson et al. (2020) Child")

dep_dict <- c(bmi = "BMI", bmi_corrected = "Corrected BMI",
              fat_ratio  = "Fat Ratio", fat_mass  = "Fat Mass", lean_mass = "Lean Mass",
              height = "Height", weight = "Weight")

save(cbbPalette, gwas_dict, file = "Data/helpers.Rdata")


# 2. Quintiles of PRS ----
get_means <- function(prs_var, dep_var, age){
  df_desc %>%
    filter(age == !!age) %>%
    select(prs = all_of(prs_var), dep_var = all_of(dep_var)) %>%
    drop_na() %>%
    mutate(prs = cut_number(prs, 5)) %>%
    lm(dep_var ~ -1 + prs, .) %>%
    tidy(conf.int = TRUE) %>%
    mutate(quintile = row_number()) %>%
    select(quintile, beta = 2, lci = 6, uci = 7)
}

df_quint <- expand_grid(prs_var = str_subset(names(df_desc), "prs"),
                        dep_var = str_subset(names(df_desc), "^(bmi|height|weight|fat|lean)"),
                        age = unique(df_desc$age)) %>%
  filter(!(dep_var %in% c("fat_ratio", "fat_mass", "lean_mass") & age != 63)) %>%
  mutate(res = pmap(list(prs_var, dep_var, age), get_means)) %>%
  unnest(res) %>%
  mutate(prs_clean = factor(gwas_dict[prs_var], gwas_dict),
         quintile = factor(quintile) %>% fct_rev(),
         dep_clean = factor(dep_dict[dep_var], dep_dict))

df_quint %>%
  filter(dep_var %in% c("fat_ratio", "fat_mass", "lean_mass")) %>%
  ggplot() +
  aes(x = quintile, y = beta, ymin = lci, ymax = uci) +
  facet_grid(prs_clean ~ dep_clean, scales = "free", switch = "y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange() +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0)) +
  labs(x = "Quintile", y = "Mean")

# 3. Kernel Density Plots ----
# BMI
plot_bmi <- function(df){
  df_kern <- df %>%
    drop_na(matches("prs"), bmi)
  
  df_stat <- df_kern %>%
    group_by(age_f) %>%
    descr(bmi) %>%
    tb() %>%
    mutate(across(c(mean, sd, skewness), round, 1),
           n = format(n.valid, big.mark = ",") %>% trimws(),
           string = glue("N  = {n}\nMean = {mean}\nSD = {sd}\nSkew = {skewness}")) %>%
    select(age_f, string)
  
  ggplot(df_kern) +
    aes(x = bmi) +
    facet_wrap(~ age_f, ncol = 3) +
    geom_density(data = select(df_kern, -age_f), aes(group = age),
                 color = "grey70", fill = "grey70", alpha = 0.4) +
    geom_density(aes(color = age_f, fill = age_f), alpha = 0.7) +
    geom_text(aes(x = Inf, y = Inf, label = string),
              data = df_stat, hjust = 1.1, vjust = 1.1) +
    guides(color = "none", fill = "none") +
    theme_bw() +
    labs(x = "BMI", y = "Density")
}

plot_bmi(df_desc)
ggsave("Images/density_obs.png",
       height = 21, width = 21, units = "cm")

df_desc %>%
  filter(sample == "cc") %>%
  plot_bmi() 
ggsave("Images/density_cc.png",
       height = 21, width = 21, units = "cm")


# PRS
plot_prs <- function(prs_var){
  df_kern <- df_desc %>%
    drop_na(matches("prs"), bmi) %>%
    select(id, age_f, age, matches("prs"), bmi) %>%
    pivot_longer(matches("prs"), names_to = "prs_var", values_to = "prs") %>%
    filter(prs_var == !!prs_var)
  
  df_stat <- df_kern %>%
    group_by(age_f, prs_var) %>%
    descr(prs) %>%
    tb() %>%
    mutate(across(c(mean, sd, skewness), round, 1),
           n = format(n.valid, big.mark = ",") %>% trimws(),
           string = glue("N  = {n}\nMean = {mean}\nSD = {sd}\nSkew = {skewness}")) %>%
    select(prs_var, age_f, string)
  
  p <- ggplot(df_kern) +
    aes(x = prs) +
    facet_wrap(~ age_f, ncol = 3) +
    geom_density(data = select(df_kern, -age_f), aes(group = age),
                 color = "grey70", fill = "grey70", alpha = 0.4) +
    geom_density(aes(color = age_f, fill = age_f), alpha = 0.7) +
    geom_text(aes(x = Inf, y = Inf, label = string),
              data = df_stat, hjust = 1.1, vjust = 1.1) +
    guides(color = "none", fill = "none") +
    theme_bw() +
    labs(x = "Polygenic Score", y = "Density")
  
  glue("Images/density_{prs_var}.png") %>%
    ggsave(plot = p, height = 21, width = 21, units = "cm")
  
  return(p)
}

names(gwas_dict) %>% map(plot_prs)


# 4. PRS, BMI, Height, Weight Correlations ----
# PRS-BMI Scatter Plot
df_scat <- df_desc %>%
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
df_desc %>%
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
       height = 21, width = 21, units = "cm")

# BMI and Height/Weight Correlations
df_desc %>%
  select(age, bmi, weight, height) %>%
  drop_na() %>%
  arrange(age) %>%
  group_split(age, .keep = FALSE) %>%
  map_dfr(~ correlate(.x, quiet = TRUE) %>%
            stretch(),
          .id = "age") %>%
  drop_na() %>%
  mutate(age = unique(df_desc$age)[as.numeric(age)] %>%
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

df_desc %>%
  distinct(id, across(matches("prs_"))) %>%
  drop_na() %>%
  select(-id) %>%
  correlate()

# 5. Attrition ----
# Distribution of BMI by PRS Missing/Observed
df_attrit <- df_desc %>%
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

distinct(df_attrit, prs_var) %>%
  expand_grid(x_var = c("bmi", "prs")) %$%
  map2(x_var, prs_var, plot_attrit)

df_desc %>%
  select(id, age, matches("prs"), bmi) %>%
  mutate(obs_bmi = ifelse(is.na(bmi), 0, 1)) %>%
  pivot_longer(matches("prs"), values_to = "prs", names_to = "prs_var") %>%
  nest(data = -c(age, prs_var)) %>%
  mutate(res = map(data, 
                   ~ lm(prs ~ obs_bmi, .x) %>%
                     tidy(conf.int = TRUE) %>%
                     filter(term == "obs_bmi") %>%
                     select(beta = estimate, lci = conf.low, uci = conf.high))) %>%
  unnest(res) %>%
  select(-data) %>%
  mutate(age = factor(age) %>% fct_rev(),
         prs_clean = factor(gwas_dict[prs_var], gwas_dict)) %>%
  ggplot() +
  aes(x = age, y = beta, ymin = lci, ymax = uci) +
  facet_wrap(~ prs_clean, nrow = 1, labeller = label_wrap_gen()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange() +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0)) +
  coord_flip() +
  labs(x = "Age", y = "Difference in PRS (Observed vs Missing BMI)")
ggsave("Images/attrit_means.png", 
       height = 21, width = 29.7, units = "cm")

rm(df_attrit)

# Difference in odds PRS observed by BMI
observed_prs <- df_desc %>% 
  select(id, matches("prs")) %>% 
  distinct() %>%
  drop_na() %>%
  pull(id)

attrit_glm <- df_desc %>%
  group_by(age) %>%
  mutate(bmi_std = wtd_scale(bmi),
         obs_prs = ifelse(id %in% observed_prs, 1, 0)) %>%
  select(id, age, obs_prs, bmi_std) %>%
  ungroup() %>%
  drop_na() %>%
  nest(data = -age) %>%
  mutate(map_dfr(data,
                 ~ glm(obs_prs ~ bmi_std, binomial, .x) %>%
                   tidy(conf.int = TRUE, exponentiate = TRUE) %>%
                   filter(term == "bmi_std") %>%
                   select(beta = estimate, lci = 6, uci = 7))) %>%
  select(-data) %>%
  mutate(age = factor(age),
         age_f = fct_rev(age))

ggplot(attrit_glm) +
  aes(x = age_f, y = beta, ymin = lci, ymax = uci) +
  geom_hline(yintercept = 1) +
  geom_pointrange() +
  scale_y_log10() +
  coord_flip() +
  theme_minimal() +
  labs(x = "Age", y = "Odds Ratio (Observed vs Missing PRS Score)")
ggsave("Images/bmi_by_prs_observed.png", 
       height = 16, width = 21, units = "cm")


# Difference in BMI by PRS Observed
attrit_lm <- df_desc %>%
  group_by(age) %>%
  mutate(bmi_std = wtd_scale(bmi),
         miss_prs = ifelse(id %in% observed_prs, 0, 1)) %>%
  select(id, age, miss_prs, bmi_std) %>%
  ungroup() %>%
  drop_na() %>%
  nest(data = -age) %>%
  mutate(map_dfr(data,
                 ~ lm(bmi_std ~ miss_prs, .x) %>%
                   tidy(conf.int = TRUE) %>%
                   filter(term == "miss_prs") %>%
                   select(beta = estimate, lci = 6, uci = 7))) %>%
  select(-data) %>%
  mutate(age = factor(age),
         age_f = fct_rev(age))

ggplot(attrit_lm) +
  aes(x = age_f, y = beta, ymin = lci, ymax = uci) +
  geom_hline(yintercept = 0) +
  geom_pointrange() +
  coord_flip() +
  theme_minimal() +
  labs(x = "Age", y = "Difference in BMI (Missing vs Observed PRS Score)")
ggsave("Images/bmi_by_prs_missing.png", 
       height = 16, width = 21, units = "cm")



# 6. PRS x Socio-Economic Position -----
df_sep <- df_desc %>%
  distinct(id, sep = father_class, across(matches("prs"))) %>%
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
  guides(color = "none", fill = "none")
ggsave("Images/prs_density_sep.png",
       height = 21, width = 29.7, units = "cm")

df_sep %>%
  nest(data = -gwas_clean) %>%
  mutate(res = map(data,
                   ~ lm(prs_value ~ - 1 + sep, .x) %>%
                     tidy(conf.int = TRUE))) %>%
  unnest(res) %>%
  select(-data) %>%
  filter(str_detect(term, "^sep")) %>%
  mutate(term = str_replace(term, "^sep", "") %>%
           ordered(levels(df_sep$sep)) %>%
           fct_rev()) %>%
  ggplot() +
  aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high,
      color = gwas_clean, shape = gwas_clean) +
  # facet_wrap(~ gwas_clean) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_pointrange(position = position_dodge(0.5)) +
  scale_shape_manual(values = 15:19) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position = "bottom") +
  coord_flip() +
  labs(x = NULL, y = "Mean Polygenic Risk Score (+95% CI)",
       color = NULL, shape = NULL) +
  guides(color = guide_legend(nrow = 2, byrow=TRUE),
         shape = guide_legend(nrow = 2, byrow=TRUE))
ggsave("Images/mean_prs_sep.png",
       height = 16, width = 21, units = "cm")

sep_means <- df_sep %>% 
  group_by(gwas_clean, sep) %>%
  summarise(prs_value = mean(prs_value), .groups = "drop")

df_sep %>%
  mutate(sep = fct_rev(sep)) %>%
  ggplot() +
  aes(x = sep, y = prs_value) +
  facet_wrap(~ gwas_clean) +
  geom_jitter(alpha = 0.3, color = cbbPalette[2]) +
  geom_boxplot(fill = NA, color = "grey50") +
  geom_violin(fill = NA, color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_point(data = sep_means, shape = 5, color = "grey50", 
             size = 3, stroke = 1) +
  theme_bw() +
  coord_flip() +
  labs(x = NULL, y = "Polygenic Risk Score")
ggsave("Images/violin_prs_sep.png",
       height = 16, width = 21, units = "cm")

rm(sep_means)


# 7. Probability of Superiority ----
set.seed(1)
df_desc %>%
  uncount(ifelse(sample == "cc", 2, 1), .id = "n") %>%
  mutate(sample = ifelse(n == 2, "obs", sample)) %>%
  select(matches("prs"), bmi, female, age, sample) %>%
  pivot_longer(matches("prs"), names_to = "prs_var", values_to = "prs") %>%
  drop_na() %>%
  group_by(age, female, prs_var, sample) %>%
  sample_n(1e4, replace = TRUE) %>%
  ungroup() %>%
  mutate(odd = ifelse(row_number() %% 2 == 0, 1, 2)) %>% 
  pivot_wider(id_cols = c(female, age, prs_var, sample),
              names_from = odd,
              names_glue = "{.value}_{odd}",
              values_from = c(bmi, prs)) %>%
  unnest(c(bmi_1, bmi_2, prs_1, prs_2)) %>%
  mutate(bmi_diff = ifelse(prs_1 > prs_2, bmi_1 - bmi_2, bmi_2 - bmi_1)) %>%
  group_by(age, prs_var, sample) %>%
  summarise(ps = sum(bmi_diff > 0)/n(),
            .groups = "drop") %>%
  mutate(age = factor(age),
         sample = ifelse(sample == "cc", "Complete Cases", "Observed Sample") %>%
           fct_rev(),
         gwas_clean = factor(gwas_dict[prs_var], gwas_dict)) %>%
  ggplot() +
  aes(x = age, y = ps, group = gwas_clean,
      color = gwas_clean, shape = gwas_clean) +
  facet_wrap(~ sample) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_line() +
  geom_point() +
  scale_shape_manual(values = 15:19) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Age", y = "Probability of Superiority",
       color = NULL,
       shape = NULL)
ggsave("Images/prob_superiority.png",
       height = 9.9, width = 21, units = "cm")

# 8. Descriptive Table ----
load("Data/height_power.Rdata")

flx <- df_desc %>%
  drop_na(matches("prs"), bmi) %>%
  select(age, bmi) %>%
  group_by(age) %>%
  descr() %>%
  tb() %>%
  left_join(height_power, by = "age") %>%
  select(age, n = n.valid, mean, sd, skewness, min, max, `power correction` = power) %>%
  rename_with(str_to_title) %>%
  rename(SD = Sd) %>%
  flextable() %>%
  border_remove() %>%
  border_inner_h(border = fp_border(color = "grey50", width = 1, style = "dashed"),
                 part = "body") %>% 
  hline_top(border = fp_border(color = "black", width = 2), part = "all") %>% 
  hline_bottom(border = fp_border(color = "black", width = 2), part = "all") %>%
  fix_border_issues(part = "all") %>% 
  align(align = "center", part = "all") %>% 
  valign(valign = "center", part = "all") %>%
  colformat_double(j = 3:7, digits = 2) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>% 
  autofit()
flx
save_as_docx(flx, path = "Tables/descriptives.docx")


# 9. PRS Correlations ----
df_desc %>%
  select(id, matches("prs")) %>%
  distinct() %>%
  select(-id) %>%
  correlate() %>%
  shave(upper = TRUE)
