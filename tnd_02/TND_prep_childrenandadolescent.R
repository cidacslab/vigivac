library(data.table)
library(tibble)
library(dplyr)
library(magrittr)
library(plotly)
library(arrow)
library(tidyverse)
library(tidylog)
library(tidytable)
library(lubridate)
library(summarytools)
library(gtsummary)


tnd=data.table::rbindlist(lapply(Sys.glob("/Users/pilarveras/Desktop/tnd19042022_gestantes_03052022.parquet"), arrow::read_parquet))

##STEP 0 - Creating variable of valid collection date (dt_coleta_valid)
##This step here was included only because we had problems with collection date in or database (2022-03-01)
##The correction was performed based on the mistakes found in the years of dt_coleta. 
##Therefore, this step, and the use of this variable is only necessary, if dt_coleta in our database 
##continues to present this problem. 

tnd = tnd %>% mutate(dt_coleta_valid=case_when(
  dt_inicio_sintomas<=dt_coleta ~ dt_coleta,
  dt_inicio_sintomas>dt_coleta  ~ dt_inicio_sintomas))



##STEP 01 - Creating variable of first case date. This step must be performed before filtering the period of study, 
##Creating the variable of first infection at this point is essential 
## so we don't lose the infections that occurred before the Study Period.

df_pos <- tnd %>% 
  filter(confirmado == "Confirmado") %>%
  select(id_vigvac, dt_coleta_valid) %>%
  arrange(dt_coleta_valid) %>%
  group_by(id_vigvac) %>%
  dplyr::mutate(counter = row_number()) %>%
  filter(counter == 1) %>%
  mutate(firstcasedate = dt_coleta_valid) %>%
  select(id_vigvac, firstcasedate)

tnd = left_join(tnd, df_pos, by = "id_vigvac")


tnd %>% count(firstcasedate,sort=T) %>% head(5)

##STEP 02 - Filtering Study Period based on collection date and date max of database 

max_date <- as.Date("2022-04-19") ##Última data do banco de vacinação 

tnd_v1 = tnd %>% filter(dt_coleta_valid>"2021-04-26" & dt_coleta_valid < max_date) ##COmeço da vacinação em criança
tnd_v2 = tnd_v1 %>% filter(assintomatico==0)
tnd_v2 = tnd_v2 %>% filter(non_symptom_info==0)

tnd_v1 %>% count(non_symptom_info,sort=T) %>% View
tnd_v1 %>% count(assintomatico,sort=T) %>% View

## quantos entram removendo o diff valid 
## diferença nos assintomaticos nos dois casos 
## inconscistencia de data (?)

##STEP 03 - Checking if there are still register that are not PCR or antigen (this step is already done, 
##in TND_prep_00- STEP 00 - SG_long and STEP 06 - SRAG_long)
#Also, check if there are no missing values in the date of first symptoms 

tnd_v2 = tnd_v2 %>% filter(tipo_teste == "RT-PCR" | tipo_teste =="Antígeno")
tnd_v2 = tnd_v2 %>% filter(!is.na(dt_inicio_sintomas))

##STEP 04 - Creating the outcome variable 
##Basically, hospitalization may occurs 4 days prior of collection date 
##up to 14 days of collection date (considering that individuals may be hospitalized with symptoms and perform the test in the hospital). 
##Death events may occur from day 0 to day 28 since date of collection. 

tnd_v3 = tnd_v2 %>%
  arrange.(dt_coleta_valid, .by=id_vigvac) %>%
  mutate.(
    diff_sample_hosp = (dt_interna - dt_coleta_valid),
    diff_sample_death = (dt_evolucao - dt_coleta_valid)
  ) %>%
  mutate.(
    hosp_event = case_when.(
      diff_sample_hosp>(-4) & diff_sample_hosp<= 14 ~ "Hospitalization" 
    ),
    death_event = case_when.(
      evolucao %in% c(2, 3) & (diff_sample_death >= 0 & diff_sample_death <= 28) ~ "Death-related"
    )
  ) %>%
  mutate.(outcome = case_when.(
    (hosp_event == "Hospitalization" | death_event == "Death-related") ~ "Hosp_Death",
    TRUE ~ "No event"
  ))



tnd_v3 %>% count(outcome,sort=T) %>% head(5)


##STEP 05 - Identifying and excluding - date inconsistencies in our database 
##Here we are excluding if: 
## - death occurred prior to the collection date 
## - hospitalization occurred after the last date of our database (max_date)
## - outcome date (dt_evolucao) occurred after max_date 
## - date of hospitalization occurred more than 14 day before first date symptoms 
## - date of collection occurred more than 14 day before first date symptoms 

tnd_v4 = tnd_v3 %>% mutate(incon_data=case_when(
  diff_sample_death<0 & evolucao %in% c(2,3)~FALSE,
  dt_interna>max_date~ FALSE,
  dt_evolucao>max_date~FALSE,
  dt_interna+14<dt_inicio_sintomas~FALSE,
  dt_coleta_valid+14<dt_inicio_sintomas~FALSE,
  TRUE~TRUE))

tnd_v4 %>% count(incon_data,sort=T) %>% head(5)
tnd_v5 = tnd_v4 %>% filter(incon_data == TRUE)


+##STEP 06 - Excluding tests representing the same infection and outcomes
  ##Here we are excluding if: 
  ## - Positive test which occurred in the last 90 days (could be the same infection)
  ## - Negative followed by positive test occurred in the last 7 days 
  ## - Negative test followed by negative test occurred in the last 14 days
  ## - For positive tests and outcome hosp_death excluded negative test occurred less then 28 days after 
  ## - For positive tests and outcome hosp_death excluded another confirmed test occurred less then 28 days after 
  
  tnd_v6 = tnd_v5 %>% group_by(id_vigvac)%>% arrange(dt_coleta_valid)%>%
  mutate(valid_test = case_when(
    confirmado == "Confirmado" & lag(confirmado == "Confirmado") & (dt_coleta_valid - lag(dt_coleta_valid)) < 90 ~ "non_valid_90",
    confirmado == "Confirmado" & lag(confirmado == "Confirmado", n = 2) & (dt_coleta_valid - lag(dt_coleta_valid, n = 2)) < 90 ~ "non_valid_90",
    confirmado == "Negativo" & lag(confirmado == "Confirmado") & (dt_coleta_valid - lag(dt_coleta_valid)) < 90 ~ "non_valid_90",
    confirmado == "Negativo" & lag(confirmado == "Confirmado", n = 2) & (dt_coleta_valid - lag(dt_coleta_valid,n=2)) < 90 ~ "non_valid_90",
    confirmado == "Negativo" & lead(confirmado == "Confirmado") & (lead(dt_coleta_valid) - dt_coleta_valid) < 7 ~ "non_valid_7",
    confirmado == "Negativo" & lead(confirmado == "Confirmado",n=2) & (lead(dt_coleta_valid,n=2) - dt_coleta_valid) < 7 ~ "non_valid_7",
    confirmado == "Negativo" & lag(confirmado == "Negativo") & (dt_coleta_valid - lag(dt_coleta_valid)) < 14 ~ "non_valid_14",
    TRUE ~ "valid"
  ))
tnd_v6 = tnd_v6 %>% ungroup()

tnd_v6 %>% count(valid_test,sort=T) %>% head(5)

tnd_v7 = tnd_v6 %>% group_by(id_vigvac)%>%
  arrange(dt_coleta_valid) %>%
  mutate(valid_outcome = (case_when(
    n() > 1 & outcome == "Hosp_Death" & lag(confirmado) == "Confirmado" & confirmado == "Negativo" & (dt_coleta_valid - lag(dt_coleta_valid)) < 28 ~ "non-valid",
    n() > 1 & outcome == "Hosp_Death" & lead(confirmado) == "Confirmado" & confirmado == "Negativo" & (lead(dt_coleta_valid) - dt_coleta_valid) < 28 ~ "non-valid",
    TRUE ~ "valid"
  )))

tnd_v7 = tnd_v7 %>% ungroup()
tnd_v7 %>% count(valid_outcome,sort=T) %>% head(5)

tnd_v8 = tnd_v7 %>% filter(valid_test=="valid")
tnd_v8 = tnd_v8 %>% filter(valid_outcome=="valid")



##STEP 07 - Excluding registers without city (cod_ibge) and sex information 

tnd_v9 = tnd_v8 %>% filter(!is.na(cod_ibge))
tnd_v10 = tnd_v9 %>% filter(!is.na(sexo) & !(sexo %in% c("I"))) 

##STEP 08 - Creating Variant of Concern (VOC) by date of prevalence  
tnd_v11 = tnd_v10 %>% mutate.(VOC=case_when.(dt_coleta_valid>="2022-01-01"~"Omicron",
                                             dt_coleta_valid<"2022-01-01"~"Delta"))
tnd_v11 %>% count(VOC,sort=T) %>% head(5)


##STEP 09 - Creating variable of vaccination status 

#FOR ADOLESCENTS
tnd_v12 = tnd_v11 %>% mutate.(status_vacinal = case_when.(
  (dt_coleta_valid - d2_data_aplicacao) > 13 ~ "full_vac",
  ((dt_coleta_valid - d1_data_aplicacao) > 13 &
     (dt_coleta_valid - d2_data_aplicacao) <= 13)|
    ((dt_coleta_valid - d1_data_aplicacao) > 13 &
       is.na(d2_data_aplicacao))~ "part_vac",
  (dt_coleta_valid - d1_data_aplicacao) >=0 &
    (dt_coleta_valid - d1_data_aplicacao) <= 6 ~ "v1_0:6",
  (dt_coleta_valid - d1_data_aplicacao) >=7 &
    (dt_coleta_valid - d1_data_aplicacao) <= 13 ~ "v1_7:13",
  (dt_coleta_valid - d1_data_aplicacao)<0 | is.na(d1_data_aplicacao) ~"uv"))

tnd_v12 %>% count(status_vacinal,sort=T) %>% head(5)


##FOR CHILDREN
tnd_v12 = tnd_v11 %>% mutate.(status_vacinal = case_when.(
  (dt_coleta_valid - d2_data_aplicacao) > 13 ~ "full_vac",
  ((dt_coleta_valid - d1_data_aplicacao) > 13 &
     (dt_coleta_valid - d2_data_aplicacao) <= 13)|
    ((dt_coleta_valid - d1_data_aplicacao) > 13 &
       is.na(d2_data_aplicacao))~ "part_vac",
  (dt_coleta_valid - d1_data_aplicacao) >=0 &
    (dt_coleta_valid - d1_data_aplicacao) <= 13 ~ "v1_0:13",
  (dt_coleta_valid - d1_data_aplicacao)<0 | is.na(d1_data_aplicacao) ~"uv"))

##For Pregnant Women

tnd_v12 = tnd_v11%>% mutate.(status_vacinal = case_when.(
  (dt_coleta_valid - d3_data_aplicacao) >= 0 ~ "v3_14",
  (dt_coleta_valid - d2_data_aplicacao) > 13 & is.na(d3_data_aplicacao) |
    (dt_coleta_valid - d2_data_aplicacao) > 13 & (dt_coleta - d3_data_aplicacao) < 0 ~ "v2_14",
  (dt_coleta_valid - d2_data_aplicacao) >= 0 & (dt_coleta - d2_data_aplicacao) <= 13 ~ "v2_0:13",
  ((dt_coleta_valid - d1_data_aplicacao) > 13 &
     (dt_coleta_valid - d2_data_aplicacao) <= 13) |
    ((dt_coleta_valid - d1_data_aplicacao) > 13 &
       is.na(d2_data_aplicacao))~ "v1_14+",
  (dt_coleta_valid - d1_data_aplicacao) >= 0 &
    (dt_coleta_valid - d1_data_aplicacao) <= 13 ~ "v1_0:13",
  (dt_coleta - d1_data_aplicacao) < 0 | is.na(d1_data_aplicacao) ~"uv"))



##STEP 10 - Creating vaccination status through length of time post first and second dose (for waning analysis)

##FOR ADOLESCENTS
tnd_v13 <- tnd_v12 %>% mutate.(
  days_vacc_1_test = as.numeric(dt_coleta_valid - d1_data_aplicacao),
  days_vacc_2_test = as.numeric(dt_coleta_valid - d2_data_aplicacao),
  status_waning = case_when.(
    days_vacc_2_test >= 0 & days_vacc_2_test <= 13 ~ "v2_0:13",
    days_vacc_2_test >= 14 & days_vacc_2_test <= 27 ~ "v2_14:27",
    days_vacc_2_test >= 28 & days_vacc_2_test <= 41 ~ "v2_28:41",
    days_vacc_2_test >= 42 & days_vacc_2_test <= 55 ~ "v2_42:55",
    days_vacc_2_test >= 56 & days_vacc_2_test <= 69 ~ "v2_56:69",
    days_vacc_2_test >= 70 & days_vacc_2_test <= 83 ~ "v2_70:83",
    days_vacc_2_test >= 84 & days_vacc_2_test <= 97 ~ "v2_84:97",
    days_vacc_2_test >= 98 ~ "v2_90+",
    days_vacc_1_test >= 0 & days_vacc_1_test <= 6 ~ "v1_0:6",
    days_vacc_1_test >= 7 & days_vacc_1_test <= 13 ~ "v1_7:13",
    days_vacc_1_test >= 14 & (days_vacc_2_test < 0 | is.na(days_vacc_2_test)) ~ "v1_p14+",
    is.na(days_vacc_1_test) | days_vacc_1_test < 0 ~ "uv"))

##FOR CHILDREN
tnd_v13 <- tnd_v12 %>% mutate.(
  days_vacc_1_test = as.numeric(dt_coleta_valid - d1_data_aplicacao),
  days_vacc_2_test = as.numeric(dt_coleta_valid - d2_data_aplicacao),
  status_waning = case_when.(
    days_vacc_2_test >= 0 & days_vacc_2_test <= 13 ~ "v2_0:13",
    days_vacc_2_test >= 14 ~ "v2_14+",
    days_vacc_1_test >= 0 & days_vacc_1_test <= 13 ~ "v1_0:13",
    days_vacc_1_test >= 14 & (days_vacc_2_test < 0 | is.na(days_vacc_2_test)) ~ "v1_p14+",
    is.na(days_vacc_1_test) | days_vacc_1_test < 0 ~ "uv"))



##In case you want to increase the lenght of time, Include: 
##days_vacc_2_test >= 98 & days_vacc_2_test <= 111 ~ "v2_14:15",
###days_vacc_2_test >= 112 & days_vacc_2_test <= 125 ~ "v2_16:17"

##STEP 11 - Excluding (for adolescents) registers with third dose register
##This step is only necessary if you don't want to check booster those yet.
##If you do want to include booster, you must treat differently your vaccination status. 

tnd_v14=tnd_v13 %>% filter(is.na(d3_nome_vacina)&is.na(d3_data_aplicacao))


##STEP 12 - Write table with all vaccines of TND participant selection 


saveRDS(tnd_v12, "/tnd_v12.rds")
saveRDS(tnd_v14, "/tnd_v14.rds")

##TND_pep02
allvac=readRDS("/tnd_v12.rds")
names(allvac)


##Including variable of epidemiological week for each date of collection 

allvac = allvac %>% mutate(SE=epiweek(dt_coleta_valid))

##Inly for childrens analysis 
allvac = allvac %>% mutate(mes=month(dt_coleta_valid))




## Change class of variables/create number of comorbididies and previous infection  
## Here we don't include pregnancy and post-part to the list of comorb. - they go separately in our analysis.  

allvac = allvac %>% mutate.(
  days_first_positive=dt_coleta_valid - firstcasedate
) %>%
  mutate_rowwise.(n_comorb = sum(c_across.(diabetes:drc))) %>%
  mutate(prev_infect = case_when(
    days_first_positive > 180 ~ ">6 mo",
    days_first_positive > 90 & days_first_positive <= 180 ~ "3-6 mo",
    days_first_positive<=90 ~ "Not previous infected",
    is.na(firstcasedate) ~ "Not previous infected"))%>% mutate.(
      n_comorb = fct_lump(factor(n_comorb), n = 2, other_level = "2+"),
      days_vacc_from_firstpos = as.numeric(d1_data_aplicacao - firstcasedate),
      raca = fct_explicit_na(factor(raca))
    )

##Only adolescents
allvac %>% count(gestante,sort=T) %>% head(5)
allvac %>% count(puerpera,sort=T) %>% head(5)


##Including variable of Brazilian Deprivation Index

ibp = fread("/Users/pilarveras/Desktop/VE-teen_children/Bancos/IBP.csv")
ibp = ibp %>% select(nome_grande_regiao,uf,codmun6,q_measure_1f_12)
allvac$cod_ibge = as.numeric(allvac$cod_ibge)
colnames(allvac)[11]  = "codmun6"
allvac = allvac %>% left_join(ibp, by="codmun6")


##Creating vs_type based in 1st dose 

allvac = allvac %>% mutate.(vs_type = paste(d1_nome_vacina, status_vacinal, sep = "_"))
allvac <- allvac %>%
  mutate.(
    vs_type = case_when.(
      vs_type == "Ad26_uv" ~ "uv",
      vs_type == "CV_uv" ~ "uv",
      vs_type == "BNT162b2_uv" ~ "uv",
      vs_type == "AZ_uv" ~ "uv",
      vs_type == "NA_uv" ~ "uv",
      TRUE ~ vs_type
    )
  )


##Children analysis CoronaVac
cv_fix = allvac %>% filter(str_detect(vs_type, "CV|uv"))
cv_fix %>% count(VOC,sort=T) %>% head(15)

bnt_fix = allvac %>% filter(str_detect(vs_type, "BNT|uv"))
bnt_fix %>% count(vs_type,sort=T) %>% head(10)



##CV - Children - Exclude tests which 1st and 2nd dose happened with less than 14 days 
cv_fix$d1_data_aplicacao=as.Date(cv_fix$d1_data_aplicacao)
cv_fix$d2_data_aplicacao=as.Date(cv_fix$d2_data_aplicacao)

cv_fix = cv_fix %>% filter(d1_data_aplicacao>="2022-01-21" | is.na(d1_data_aplicacao)|is.na(d2_data_aplicacao))



##BNT - Adolescent -  Exclude tests which 1st and 2nd dose happened with less than 14 days 
bnt_fix$d1_data_aplicacao=as.Date(bnt_fix$d1_data_aplicacao)
bnt_fix$d2_data_aplicacao=as.Date(bnt_fix$d2_data_aplicacao)

bnt_fix = bnt_fix %>% mutate.(diff_vac = d2_data_aplicacao - d1_data_aplicacao) %>% 
  filter(diff_vac>=14 | is.na(d1_data_aplicacao)|is.na(d2_data_aplicacao))


saveRDS(cv_fix, "/Users/pilarveras/Desktop/TND_CV_19042022_6a11_25042022.rds")

saveRDS(bnt_fix, "/Users/pilarveras/Desktop/TND_BNT_19042022_5a11_25042022.rds")