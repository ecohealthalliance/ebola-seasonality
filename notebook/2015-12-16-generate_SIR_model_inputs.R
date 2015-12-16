# Load in and reshape the data

source('../R/process_sheet.R') #load the sheet-processing function
raw_xl = '../data/P giganteus Filo serology_ARW_NR_JE_ARW_12.11.xls'
data_sheets = excel_sheets(raw_xl) %>%`[`(stri_detect_regex(., "\\d{4}")) %>% `[`(stri_detect_regex(., "Faridpur"))
data_list = lapply(data_sheets, process_sheet, excel_file = raw_xl)
bat_serol = bind_rows(data_list) %>%
  gather("antibody", "luminex_val", matches("(sg|sf|vn|gp|_scale)$")) %>%
  mutate(tmp = ifelse(stri_detect_fixed(antibody, '_scale'), 'sero', 'luminex'),
         antibody = stri_replace_first_fixed(antibody, '_scale', '')) %>%
  spread(tmp, luminex_val) %>%
  mutate(luminex = as.numeric(luminex))

#  Use only ebola data, calculate seroprevalence

filtered_averaged = bat_serol %>%
  filter(antibody == "ebolagp") %>%
  group_by(collection_month, year) %>%
  mutate(collection_date = round_date(mean(date_corrected), "day")) %>%
  group_by() %>%
  #  mutate(age_class = ifelse(age_class == "preweaned", "juvenile", age_class)) %>%
  filter(age_class != "preweaned") %>%
  filter(!is.na(age_class)) %>%
  group_by(collection_date, age_class) %>%
  summarize(count = sum(!is.na(sero)),
            sero_pos = sum((sero != "neg") & !is.na(sero)),
            seroprevalence = sero_pos/count,
            sero_err = sqrt((seroprevalence*(1-seroprevalence))/count))

# Reshape to look like `data/niVFarlong`

old_niv_data = read.csv('../data/niVFarLong.csv')

new_eb_data = filtered_averaged %>%
  select(-sero_err) %>%
  gather("var", "val", seroprevalence, count, sero_pos) %>%
  mutate(var = paste(age_class, var, sep=".")) %>%
  select(-age_class) %>%
  spread(var, val) %>%
  mutate(Rweek = round((collection_date - ymd('2006-01-01'))/dweeks(1))) %>%
  filter(Rweek != 118) %>%
  mutate(Date = floor_date(collection_date, 'month')) %>%
  mutate(N = old_niv_data$N) %>%
  mutate(Fjuv = juvenile.count / (juvenile.count + adult.count),
         Njuv = N * Fjuv, Nadult = N * (1-Fjuv),
         Juv.tot = juvenile.count, Adult.tot = adult.count,
         Adult.pos = adult.sero_pos,
         Juv.pos = juvenile.sero_pos,
         Adult.neg = Adult.tot - Adult.pos,
         Juv.neg = Juv.tot - Juv.pos) %>%
  mutate_each(funs(as.integer), matches("\\.\\w{3}")) %>%
  select(-adult.seroprevalence, -juvenile.seroprevalence, -adult.count, -juvenile.count, -adult.sero_pos, -juvenile.sero_pos)


