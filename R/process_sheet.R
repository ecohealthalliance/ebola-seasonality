#' Import and process one data sheet of bat serology data from the Excel files
#' @export
process_sheet = function(sheet_name, excel_file) {
  require(readxl)
  require(magrittr)
  require(dplyr)
  require(tidyr)
  require(lubridate)
  require(stringi)
  require(assertr)

  capture.output(sheet <- read_excel(excel_file, sheet_name, col_names = TRUE, na = "999")[1:200,1:60], file = '/dev/null')
  sheet = sheet[, !is.na(names(sheet))]

  names(sheet) = names(sheet) %>%
    stri_trim_both() %>%
    stri_trans_tolower() %>%
    stri_replace_all_fixed(" ", "") %>%
    stri_replace_all_fixed("-", "")

  sheet = sheet %>%
    filter(!is.na(id)) %>%
    `[`(, !is.na(names(.)))



  sheet_data = sheet_name %>% stri_split_fixed(" ") %>% unlist
  sheet$region = sheet_data[1]
  sheet$collection_month = sheet_data[2]
  sheet$year = as.numeric(sheet_data[3])

  field_names = c("hevsg", "hevsg_scale", "nivsg", "nivsg_scale", "nivsf",
                  "nivsf_scale", "hevn", "hevn_scale", "nivn", "nivn_scale",
                  "cedvsg", "cedvsg_scale", "ebolagp", "ebolagp_scale", "marvgp",
                  "marvgp_scale", "menvn", "menvn_scale", "nbvn", "nbvn_scale",
                  "bcs_excellent", "microchip_id", "feces_rectal_lysis", "feces_rectal_vtm")

  for(aname in field_names) {
    if(!(aname %in% names(sheet)))
      sheet[[aname]] <- NA
  }

  sheet_processed = sheet %>%
    mutate(id = as.character(id)) %>%
    verify(stri_detect_regex(id, "(^\\d|^ENB)")) %>%
    #   verify(length(ID) == length(unique(ID))) %>%
    mutate(date_corrected = as.POSIXct(date_corrected)) %>%
    select(-date_sampled) %>%
    verify(year(date_corrected) == year) %>%
    verify((month(date_corrected, label = TRUE) == collection_month) |
             (month(date_corrected %m-% months(1), label = TRUE) == collection_month)) %>%
    mutate(feces_rectal_lysis = ifelse(stri_detect_regex(feces_rectal_lysis, "[Nn][Oo]"), 0, feces_rectal_lysis)) %>%
    mutate(feces_rectal_vtm = ifelse(stri_detect_regex(feces_rectal_vtm, "[Nn][Oo]"), 0, feces_rectal_vtm)) %>%
    verify((is.na(sex_male) != is.na(sex_female)) | (is.na(sex_male) & is.na(sex_female))) %>%
    mutate(location_id = as.character(location_id)) %>%
    mutate(sex = ifelse(is.na(sex_male) & !is.na(sex_female), "F", ifelse(is.na(sex_male) & is.na(sex_female), NA, "M"))) %>%
    select(-sex_female, -sex_male) %>%
    verify((is.na(age_preweaned) + is.na(age_juvenile) + is.na(age_adult)) >= 2) %>%
    mutate(age_class = ifelse(!is.na(age_preweaned), "preweaned", ifelse(!is.na(age_juvenile), "juvenile", "adult"))) %>%
    mutate(age_class = ifelse(is.na(age_preweaned) + is.na(age_juvenile) + is.na(age_adult) == 3, NA, age_class)) %>%
    select(-age_preweaned, -age_juvenile, -age_adult) %>%
    assert(within_bounds(50, 250), forearm_mm) %>%
    assert(within_bounds(10, 100), head_mm) %>%
    assert(within_bounds(10, 1500), weight_g) %>%
    #  verify(!(!is.na(pregnant_yes) & !is.na(pregnant_no))) %>%
    mutate(pregnant = ifelse(!is.na(pregnant_yes), TRUE, FALSE)) %>%
    # select(-pregnant_yes, -pregnant_no) %>%
    verify(!(!is.na(lactating_yes) & !is.na(lactating_no))) %>%
    mutate(lactating = ifelse(!is.na(lactating_yes), TRUE, FALSE)) %>%
    select(-lactating_yes, -lactating_no) %>%
    verify(!(!is.na(pup_yes) & !is.na(pup_no))) %>%
    mutate(pup = ifelse(!is.na(pup_yes), TRUE, FALSE)) %>%
    select(-pup_yes, -pup_no) %>%
    #    verify(!((sex == "M") & (pregnant | lactating | pup))) %>%
    verify((is.na(bcs_poor) + is.na(bcs_fair) + is.na(bcs_good) + is.na(bcs_excellent)) == 3) %>%
    mutate(body_condition = ifelse(!is.na(bcs_poor), "poor", ifelse(!is.na(bcs_fair), "fair", ifelse(!is.na(bcs_good), "good", ifelse(!is.na(bcs_excellent), "excellent", NA))))) %>%
    select(-bcs_excellent, -bcs_good, -bcs_fair, bcs_poor) %>%
    assert(in_set(1, 3, 6, 22, 36, 40, 43, 56), ends_with("_scale")) %>%
    mutate_each(funs(ifelse(. == 1, "neg", ifelse(. %in% c(36, 43), "low", ifelse(. %in% c(22,3, 6, 40), "high", NA)))), ends_with("_scale")) %>%
    verify((!stri_detect_regex(microchip_id, '[a-zA-Z]') | is.na(microchip_id))) %>%
    mutate_each(funs(as.numeric), matches("(sg|sf|vn|gp)$")) %>%
    assert(within_bounds(1, 50000), matches("(sg|sf|vn|gp)$")) %>%
    mutate(microchip_id = stri_replace_all_regex(as.character(microchip_id), "\\.0+", '')) %>%
    mutate(microchip_id = stri_replace_all_regex(as.character(microchip_id), "[^\\d]", '')) %>%
    verify((nchar(microchip_id) %in% 8:9) | is.na(microchip_id)) %>%
    mutate(microchip_id = ifelse(nchar(microchip_id) == 8, paste0("0", microchip_id), microchip_id)) %>%
    mutate(microchip_id = paste(substr(microchip_id, 1, 3), substr(microchip_id, 4, 6), substr(microchip_id, 7, 9), sep="_"))
  return(sheet_processed)
}