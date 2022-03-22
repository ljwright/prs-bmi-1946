get_sep_r2 <- function(spec_id, sep_var){
  
  df_mod <- get_df(spec_id) %>%
    rename(sep_var = all_of(!!sep_var)) %>%
    select(id, dep_var, prs, sep_var, female) %>%
    drop_na()
  
  get_boot <- function(boot){
    set.seed(boot)
    
    df_boot <- sample_frac(df_mod, replace = TRUE)
    
    get_r2 <- function(covars){
      get_form(spec_id, covars) %>%
        as.formula() %>%
        lm(df_boot) %>%
        broom::glance(mod) %>%
        pull(1)
    }
    
    get_r2(c("prs", "sep_var")) - get_r2(c("prs"))
  }
  
  map_dbl(1:500, get_boot) %>%
    get_ci()
}

res_sep_r2 <- sep_specs %>%
  filter(sex == "all", 
         prs_var == "prs_k", 
         dep_var == "bmi_raw",
         str_detect(sep_var, "ridit")) %>%
  get_furrr2(sep_var, get_sep_r2)

res_sep_r2 %>%
  mutate(age_f = factor(age)) %>%
  ggplot() +
  aes(x = age_f, y = beta, ymin = lci, ymax = uci, group = sep_var) +
  facet_grid(~ sep_var, scales = "free", switch = "y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines")) +
  labs(x = "Age", y = NULL)


df_reg %>% 
  drop_na(prs_k) %>%
  nest(data = -age) %>%
  mutate(r2 = map_dbl(data, ~ lm(bmi_raw ~ female, .) %>%
                     glance() %>% pull(1))) %>%
  arrange(age) %>%
  select(-data)





df_attrit <- clean_res(res_attrit_2) %>%
  mutate(age_from = factor(age_from) %>%
           fct_recode("Observed" = "0", "Complete Cases" = "2")) %>%
  filter(age_from != "6", age_from != "69")

ggplot(df_attrit) +
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
  labs(x = "Age", y = "Incremental R2")
