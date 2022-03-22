library(gamlss)

rm(list = ls())

# 1. Load Data ----
load("Data/df_long.Rdata")

df_reg <- df_long %>%
  group_by(age) %>%
  mutate(across(c(bmi, bmi_corrected),
                list(std = wtd_scale,
                     ln = log,
                     rank = ~ percent_rank(.x)*100))) %>%
  ungroup() %>%
  rename(bmi_raw = bmi, bmi_corrected_raw = bmi_corrected) %>%
  drop_na(matches("^prs")) %>%
  dplyr::select(age, bmi_raw, female, prs_k) %>%
  drop_na()

rm(df_long)

new_data <- expand_grid(female = c(0, 1),
                        prs_k = -2:2)


get_gamlss <- function(age){
  
}




res <- tibble()
for (age in unique(df_reg$age)){
  df_x <- df_reg %>%
    filter(age == !!age)
  
  mod <- gamlss(bmi_raw ~ female + prs_k,
                sigma.formula = ~ female + prs_k,
                nu.formula = ~ female + prs_k,
                family = BCCGo,
                data = df_x,
                trace = FALSE)
  
  pred <- predictAll(mod, newdata = new_data, type = "link") %>%
    as_tibble() %>%
    bind_cols(new_data, .)
  
  curr_res <- tibble(age = age, mod = list(mod), pred = list(pred))
  
  res <- bind_rows(res, curr_res)
}
rm(df_x, mod, pred, curr_res)

res %>%
  unnest(pred) %>%
  dplyr::select(-mod) %>%
  uncount(99, .id = "centile") %>%
  mutate(tau = centile/100,
         bmi = exp(mu) * (1 + nu * exp(sigma) * qnorm(tau))^(1/nu)) %>%
  filter(prs_k %in% c(-1, 1)) %>%
  mutate(prs_k = ifelse(prs_k == -1, "low", "high")) %>%
  pivot_wider(id_cols = c(age, female, centile), names_from = prs_k, values_from = bmi) %>%
  mutate(diff = high - low) %>%
  ggplot() +
  aes(x = centile, y = diff, color = factor(female)) +
  facet_wrap(~ age) +
  geom_hline(yintercept = 0) +
  geom_line()

