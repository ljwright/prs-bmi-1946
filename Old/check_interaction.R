sep_specs <- mod_specs %>%
  expand_grid(sep_var = c("medu", "medu_ridit", "sep_ridit", 
                          "sep_f", "sep_f4", "sep_f5"))

df_reg <- df_reg %>%
  mutate(sep_n = as.numeric(sep),
         sep_f = sep,
         sep_f4 = ifelse(between(sep_n, 1, 3), "Upper", "Lower") %>%
           factor(c("Upper", "Lower")),
         sep_f5 = ifelse(between(sep_n, 1, 4), "Upper", "Lower") %>%
           factor(c("Upper", "Lower")))



get_mult <- function(spec_id, sep_var){
  
  df_mod <- get_df(spec_id) %>%
    rename(sep_var = all_of(!!sep_var))
  
  get_form(spec_id, "prs*sep_var") %>%
    as.formula() %>%
    lm(df_mod) %>%
    broom::tidy(conf.int = TRUE) %>%
    filter(str_detect(term, "prs\\:")) %>%
    select(term, beta = 2, p = 5, lci = 6, uci = 7)
}         

res_int <- sep_specs %>%
  filter(sex == "all") %>%
  mutate(res = map2(spec_id, sep_var, get_mult)) %>%
  unnest(res) %>%
  mutate(term = str_replace(term, "prs\\:sep_var", ""),
         term = ifelse(term == "", "Interaction", term),
         age_f = factor(age))

res_int %>%
  filter(prs_var == "prs_v",
         str_detect(dep_var, "bmi|height|weight")) %>%
  filter(sep_var != "sep_f") %>%
  ggplot() +
  aes(x = age_f, y = beta, ymin = lci, ymax = uci,
      group = term, color = term, fill = term) +
  facet_wrap(sep_var ~ dep_var, scales = "free_y") +
  geom_hline(yintercept = 0) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line() +
  geom_point()

res_int %>%
  filter(prs_var == "prs_k",
         str_detect(dep_var, "bmi|height|weight")) %>%
  filter(sep_var == "sep_f") %>%
  ggplot() +
  aes(x = age_f, y = beta, ymin = lci, ymax = uci,
      group = term, color = term, fill = term) +
  facet_grid(dep_var ~ term, scales = "free_y") +
  geom_hline(yintercept = 0) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line() +
  geom_point()

res_int %>%
  filter(prs_var == "prs_k",
         !str_detect(dep_var, "bmi|height|weight")) %>%
  ggplot() +
  aes(x = term, y = beta, ymin = lci, ymax = uci,
      group = term, color = term, fill = term) +
  facet_grid(dep_var ~ sep_var, scales = "free_x") +
  geom_hline(yintercept = 0) +
  geom_pointrange() + 
  coord_flip()

get_sep <- function(spec_id, sep_var){
  df_mod <- get_df(spec_id) %>%
    rename(sep_var = all_of(!!sep_var))
  
  get_form(spec_id, "sep_var") %>%
    as.formula() %>%
    lm(df_mod) %>%
    broom::tidy(conf.int = TRUE) %>%
    filter(str_detect(term, "sep_var")) %>%
    select(term, beta = 2, p = 5, lci = 6, uci = 7)
}

res_sep <- sep_specs %>%
  filter(sex == "all", prs_var == "prs_k",
         dep_var == "bmi_raw") %>%
  mutate(res = map2(spec_id, sep_var, get_sep)) %>%
  unnest(res) %>%
  mutate(term = str_replace(term, "sep_var", ""),
         term = ifelse(term == "", "SEP", term),
         age_f = factor(age))

res_sep %>%
  ggplot() +
  aes(x = age_f, y = beta, ymin = lci, ymax = uci,
      group = term, color = term, fill = term) +
  facet_grid(term ~ sep_var, scales = "free_y") +
  geom_hline(yintercept = 0) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line() +
  geom_point()

