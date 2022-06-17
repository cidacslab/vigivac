##############################################################################
# Name of file: 02_Modelling TND
# Original author(s): Thiago Cerqueira
# Original date: 26 Feb 22
# Latest update author (if not using version control) - thiago.c.silva@fiocruz.br
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 4.1.2
# Description of content: Modelling
# Approximate run time: Unknown
##############################################################################
source("utility_functions.R")
pacman::p_load(tidyverse, lubridate, tidytable, tidylog, mgcv)
capitais <- read_csv("Capitais.csv", col_names = c("city", "cod"))
cod_ibge <- read_csv("BDI_Municipalities-Level_Short.csv")
cod_ibge <- cod_ibge %>% select(nome_grande_regiao, nome_da_uf, codmun6, q_measure_1f_12)
cod_ibge$codmun6 <- as.character(cod_ibge$codmun6)
tnd_db <- data.table::rbindlist(lapply(Sys.glob("/dados/Analysis/Thiago/Scripts/processing_data/tnd170422_no_diff_vac/part-*.parquet"), arrow::read_parquet))
tnd_db <- tnd_db %>% filter(dt_coleta_valid>"2021-12-31")
tnd_db <- tnd_db %>% filter(valid_test == "valid", valid_outcome == "valid")



z_pos <- tnd_db %>% filter(confirmado=="Confirmado")
z_neg <- tnd_db %>% filter(confirmado=="Negativo")

z_sel <- z_neg$id_vigvac %in% z_pos$id_vigvac
table(z_sel)

#remove any patients in the negative dataset that also have a positive test
z_neg <- z_neg %>% filter(!z_sel) 


# select only the first negative and first positive
z_neg <- z_neg %>%arrange.(id_vigvac,dt_coleta_valid) %>% group_by(id_vigvac) |>  filter(row_number()==1)


z_pos <- z_pos %>%arrange.(id_vigvac,dt_coleta_valid) %>%   group_by(id_vigvac) |>  filter(row_number()==1)

tnd_db <- bind_rows(z_neg,z_pos)
#

# create new variables
tnd_db <- tnd_db %>%
  mutate.(
    prev_infec = case_when.(
      dt_coleta_valid - firstcasedate > 180 ~ ">6 mo",
      dt_coleta_valid - firstcasedate > 90 ~ "3-6 mo",
      dt_coleta_valid - firstcasedate <= 90 ~ "Not",
      is.na(firstcasedate) ~ "Not"
    ),
    prev_infec = fct_relevel(
      factor(prev_infec), "Not",
      "3-6 mo",
      ">6 mo"
    )
  ) %>%
  mutate_rowwise.(n_comorb = sum(c_across.(diabetes:drc))) %>%
  mutate.(
    sexo = as.factor(sexo),
    confirmado = as.factor(confirmado),
    cod_ibge = as.factor(cod_ibge)
  ) %>%
  mutate.(
    uf = factor(str_extract(cod_ibge, "^.{2}")),
    raca = fct_recode(factor(raca),
                      White = "1",
                      Black = "2",
                      Asian = "3",
                      Mixed = "4",
                      Indigenous = "5"
    )
  ) %>%
  mutate.(
    n_comorb_cat = fct_lump(factor(n_comorb), n = 3, other_level = "3+"),
    raca = fct_explicit_na(factor(raca))
  ) %>%
  mutate.(
    capital_not = case_when.(
      cod_ibge %in% all_of(capitais$cod) ~ "Capital",
      TRUE ~ "Not"
    )
  ) %>%
  left_join(cod_ibge, by = c("cod_ibge" = "codmun6")) %>%
  mutate(q_measure_1f_12 = factor(if_else(is.na(q_measure_1f_12) & uf == 53, 2, q_measure_1f_12)))

# Create vaccine categories
tnd_db <- tnd_db %>% mutate.(
  days_vacc_1_test = as.numeric(dt_coleta_valid - d1_data_aplicacao),
  days_vacc_2_test = as.numeric(dt_coleta_valid - d2_data_aplicacao),
  days_vacc_3_test = as.numeric(dt_coleta_valid - d3_data_aplicacao),
  days_vacc_4_test = as.numeric(dt_coleta_valid - d4_data_aplicacao),
  vacc_status = case_when.(
    days_vacc_4_test >= 0 ~ "v4",
    days_vacc_3_test >= 0 & days_vacc_3_test <= 13 ~ "v3_0:13",
    days_vacc_3_test >= 14 & days_vacc_3_test <= 30 ~ "v3_14:30",
    days_vacc_3_test >= 31 & days_vacc_3_test <= 60 ~ "v3_31:60",
    days_vacc_3_test >= 61 & days_vacc_3_test <= 90 ~ "v3_61:90",
    days_vacc_3_test >= 91 & days_vacc_3_test <= 120 ~ "v3_91:120",
    days_vacc_3_test >= 121 ~ "v3_121",
    days_vacc_2_test >= 0 & days_vacc_2_test <= 13 ~ "v2_0:13",
    days_vacc_2_test >= 14 & days_vacc_2_test <= 180 ~ "v2_14:180",
    days_vacc_2_test >= 181 ~ "v2_181",
    days_vacc_1_test >= 0 & days_vacc_1_test <= 13 ~ "v1_0:1",
    days_vacc_1_test >= 14 & (days_vacc_2_test < 0 | is.na(days_vacc_2_test)) ~ "v1_2",
    is.na(days_vacc_1_test) | days_vacc_1_test < 0 ~ "uv"
  )
)

#Filter out 4 doses
tnd_db <- tnd_db |> filter(vacc_status!="v4")


# Filter out inconsistencies 
tnd_db <- tnd_db |> filter(case_when(
  str_detect(vacc_status,"v2|v3")~ no_vacc_incon& different_vaccines,
  vacc_status!="uv" ~ no_vacc_incon,
  vacc_status=="uv" ~ TRUE
))

#### People vaccinated after test- as unvaccinated
# Add vaccine type

tnd_db <- tnd_db %>%
  mutate.(vs_type = paste(d1_nome_vacina, vacc_status, sep = "_"))

tnd_db <- tnd_db %>%
  mutate.(
    vs_type = case_when.(
      vs_type == "Ad26_uv" ~ "uv",
      vs_type == "CV_uv" ~ "uv",
      vs_type == "BNT162b2_uv" ~ "uv",
      vs_type == "AZ_uv" ~ "uv",
      vs_type == "NA_uv" ~ "uv",
      TRUE ~ vs_type
    )
  ) %>%
  mutate(vs_type = fct_relevel(factor(vs_type), "uv")) %>%
  mutate(vs_type2 = case_when(
    str_detect(vs_type, "v3_") & d3_nome_vacina == "BNT162b2" ~ paste(vs_type, d3_nome_vacina, sep = "_"),
    str_detect(vs_type, "v3_") & d3_nome_vacina == "CV" ~ paste(vs_type, d3_nome_vacina, sep = "_"),
    str_detect(vs_type, "v3_") & d3_nome_vacina %in% c("Ad26","AZ") ~ paste(vs_type, "Other", sep = "_"), # combine other booster in only one category
    TRUE ~ as.character(vs_type)
  ), ) %>%
  droplevels()

# create symptomatic and severe outcomes 
tnd_db <- tnd_db %>%
  mutate(
    vs_type2 = fct_relevel(factor(vs_type2), "uv"),
    temporal_trend = as.numeric(dt_coleta_valid - as.Date("2021-12-31")),
    outcome_confirmado = case_when(
      outcome == "Hosp_Death" & confirmado == "Confirmado" ~ "Hosp_Death_Positive",
      outcome != "Hosp_Death" & confirmado == "Confirmado" ~ "Outpatient_Positive",
      confirmado == "Negativo" ~ "Negative"
    ),
    age_cut = cut(idade, breaks = c(17, seq(24, 90, 5), +Inf)),
    age_group = cut(idade, breaks = c(17, 59, 79, +Inf)),
    week_dummy = factor(paste(lubridate::epiyear(dt_coleta_valid), lubridate::epiweek(dt_coleta_valid), sep = "_")),
    across(where(is.character), as.factor)
  ) %>%
  droplevels()

# create time between doses
tnd_db <- tnd_db %>% mutate(
  diff1_2 = case_when(
    str_detect(vs_type2, "v3|v2") & d2_data_aplicacao - d1_data_aplicacao > 41 ~ ">=6wk",
    str_detect(vs_type2, "v3|v2") & d2_data_aplicacao - d1_data_aplicacao > 20 ~ "3-5wk",
    str_detect(vs_type2, "v3|v2") & d2_data_aplicacao - d1_data_aplicacao < 21 ~ "<3wk",
  ),
  diff1_2=fct_relevel(factor(diff1_2),
                      "<3wk",
                      "3-5wk",
                      ">=6wk"),
  diff2_3 = case_when(
    str_detect(vs_type2, "v3") & d3_data_aplicacao - d2_data_aplicacao > 179 ~ "180d+",
    str_detect(vs_type2, "v3") & d3_data_aplicacao - d2_data_aplicacao > 149  ~ "150-179",
    str_detect(vs_type2, "v3") & d3_data_aplicacao - d2_data_aplicacao > 114 ~ "115-149",
    str_detect(vs_type2, "v3") & d3_data_aplicacao - d2_data_aplicacao < 115 ~ "<115d",
    str_detect(vs_type2, "v2") ~ "Only 2 doses",
    str_detect(vs_type2, "v1") ~ "Only 1 dose",
    str_detect(vs_type2, "uv") ~ "Unvax",
  ),
  diff2_3=fct_relevel(factor(diff2_3),
                      "<115d",
                      "115-149",
                      "150-179",
                      "180d+",
                      "Only 2 doses")
)


# filter out different vaccines
tnd_cv <- tnd_db %>%
  filter(str_detect(vs_type2, "CV_v|uv")) |> 
  filter(str_detect(vs_type2, "Other$", negate = T)) %>%
  mutate(
    vs_type2 = fct_relevel(
      vs_type2,
      "uv",
      "CV_v1_0:1",
      "CV_v1_2",
      "CV_v2_0:13",
      "CV_v2_14:180",
      "CV_v2_181",
      "CV_v3_0:13_BNT162b2",
      "CV_v3_14:30_BNT162b2",
      "CV_v3_31:60_BNT162b2",
      "CV_v3_61:90_BNT162b2",
      "CV_v3_91:120_BNT162b2",
      "CV_v3_121_BNT162b2",
      "CV_v3_0:13_CV",
      "CV_v3_14:30_CV",
      "CV_v3_31:60_CV",
      "CV_v3_61:90_CV",
      "CV_v3_91:120_CV",
      "CV_v3_121_CV",
    )
  ) %>% filter(diff2_3 != "<115d") |> 
  filter(assintomatico==0) %>%
  droplevels()

# Remove coronavac boosters
tnd_cv_final <- tnd_cv %>% filter(str_detect(vs_type2,"CV$", negate = T)) %>% droplevels()
saveRDS(tnd_cv_final,"tnd_cv_1804.RDS")
