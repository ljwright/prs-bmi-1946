do_file <- read_lines("Data/STATA_Script_dbannZZdhxddp_20220216-115754.do") %>%
  tibble(line = .)

value_labels <- do_file %>%
  filter(str_detect(line, "^label define")) %>%
  mutate(line = str_replace(line, "label define ", "") %>%
           str_replace(", modify", "") %>%
           str_replace_all('\\"', "XX")) %>%
  separate(line, c("line", "label"), sep = " XX", extra = "drop") %>%
  mutate(label = str_replace(label, "XX", "")) %>%
  separate(line, c("var", "value"), sep = " ") %>%
  mutate(value = str_replace_all(value, "\\'", "") %>%
           as.integer())

variable_labels <- do_file %>%
  filter(str_detect(line, "^label variable")) %>%
  mutate(line = str_replace(line, "label variable ", "") %>%
           str_replace_all('\\"', "XX")) %>%
  separate(line, c("variable", "label"), sep = " XX", extra = "drop") %>%
  mutate(label = str_replace(label, "XX", ""))
