##############################################################################
# Name of file: 01.TND_Prep_data_test
# Original author(s): Thiago Cerqueira
# Original date: 11 Nov 21
# Latest update author (if not using version control) - thiago.c.silva@fiocruz.br
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 4.1.1
# Description of content: Prepares testing data
# Approximate run time: Unknown
##############################################################################
library(tidyverse)
library(tidylog)
library(tidytable)
library(tictoc)

setwd("/dados/Analysis/Thiago/Scripts/CoronaVac-Waning/TND")
tic()
lab.data <- read_rds("/dados/Analysis/Thiago/Exploratory/tnd0812.RDS")
#number of individuals 24495884
# Fill date of sample collected with notification date, if sample collected date is missing and the sample is positive  
lab.data <- lab.data %>% mutate(dt_coleta=case_when(
  is.na(dt_coleta) & confirmado=="Confirmado"~dt_notificacao,
  TRUE~dt_coleta
))
#mutate: changed 363,857 values (1%) of 'dt_coleta' (363857 fewer NA)
# Filter inconsistencies and notifications before covid-19 pandemic 
lab.data <- lab.data %>% filter(dt_coleta>"2020-02-20")
#filter: removed 12,857 rows (<1%), 30,946,602 rows remaining
max.date <- as.Date("2021-12-08")

#impute sample collected date in sympton onset, when asymptomatic
lab.data <- lab.data %>% mutate(dt_inicio_sintomas=coalesce(dt_inicio_sintomas,dt_coleta))

#mutate: changed 363,320 values (1%) of 'dt_coleta' (363320 fewer NA)

# Remove individuals <18
lab.data <- lab.data %>% filter(idade>17) 
#filter: removed 3,046,118 rows (10%), 27,900,484 rows remaining

# Get the first positive test date
df_pos <-  lab.data %>%
  filter.(confirmado=="Confirmado") %>%
  select.(id_vigvac,dt_coleta) %>%
  arrange.(id_vigvac,dt_coleta,.by=id_vigvac) %>%
  mutate.(counter = row_number(id),.by = id_vigvac) %>%
  filter.(counter==1,.by = id_vigvac) %>%
  mutate.(firstcasedate = dt_coleta) %>% 
  select.(id_vigvac,firstcasedate)

# Join first positive test date
lab.data <-  left_join(lab.data,df_pos,by="id_vigvac")
# left_join: added one column (firstcasedate)
# > rows only in x   12,740,035
# > rows only in y  (         0)
# > matched rows     15,160,449
# >                 ============
#   > rows total       27,900,484

# Remove cases before vaccine are available
lab.data <- lab.data %>% filter(dt_inicio_sintomas>"2021-01-17")
#filter: removed 8,850,298 rows (32%), 19,050,186 rows remaining

lab.data <- lab.data %>% filter(!duplicated(paste0(id_vigvac,dt_coleta,confirmado)))
#filter: removed 571,682 rows (3%), 18,478,504 rows remaining

# Check if there is a positive less than 7 days apart from a negative
lab.data <- lab.data %>%
  arrange.(id_vigvac,dt_coleta,confirmado,.by=id_vigvac) %>%
  mutate.(NegFollowedByPos = as.numeric(
    confirmado=="Negativo" &
      leads.(confirmado,default="Negativo")=="Confirmado" &
      as.numeric(leads.(dt_coleta,default=max.date+8)-dt_coleta)<=7),.by = id_vigvac)

#remove false negative ones (negative tests followed by positive) - 
lab.data <- lab.data %>% filter(NegFollowedByPos==0)
#filter: removed 95,320 rows (1%), 18,383,184 rows remaining

#check numbers of cases/controls
table(lab.data$confirmado, useNA = "always")


# Confirmado   Negativo       <NA> 
#   8396245    9986939          0 


# Create variable of date difference across tests
lab.data <- lab.data %>% 
arrange.(id_vigvac,dt_coleta) %>%
  mutate.(DaysDiffTest = dt_coleta - lags.(dt_coleta),.by = id_vigvac) 


# Fill the records of e-SUS with information from SIVEP
lab.data <- lab.data  %>% arrange.(id_vigvac,dt_interna_dt) %>% fill.(dt_interna_dt,evolucao,dt_evoluca, .direction = "down",.by = id_vigvac)

#remove consecutive tests of people with same result At least 14 days of difference between negatives and 90 for positive

lab.data <- lab.data %>%group_by(id_vigvac) %>%  filter(case_when(n()>1 & all(confirmado=="Negativo")~DaysDiffTest>14 | is.na(DaysDiffTest),
                                  TRUE~TRUE))%>%  filter(case_when(n()>1 & all(confirmado=="Confirmado")~DaysDiffTest>90 | is.na(DaysDiffTest),
                                                                 TRUE~TRUE))
# group_by: one grouping variable (id_vigvac)
# filter (grouped): removed 292,514 rows (2%), 18,090,670 rows remaining
# filter (grouped): removed 333,180 rows (2%), 17,757,490 rows remaining

lab.data <- lab.data %>% 
  arrange.(id_vigvac,dt_coleta) %>%
  mutate.(DaysDiffTest = dt_coleta - lags.(dt_coleta),.by = id_vigvac) 

# Deal with people with positive and negative. Remove positive followed by positive less than 90 days - 

lab.data <- lab.data %>%group_by(id_vigvac) %>%arrange(id_vigvac,dt_coleta) %>%  
  filter(case_when(n()>1 & lag(confirmado)=="Confirmado" & confirmado=="Confirmado"~DaysDiffTest>90,
                   n()>2 & lag(confirmado,n=2)=="Confirmado" & confirmado=="Confirmado"~dt_coleta -lag(dt_coleta,n=2)>90,
                   n()>1 & confirmado=="Negativo" & lead(confirmado)=="Confirmado"~lead(dt_coleta)-dt_coleta>7 | is.na(DaysDiffTest),
                   TRUE~TRUE)) 


# group_by: one grouping variable (id_vigvac)
# filter (grouped): removed 79,180 rows (<1%), 17,678,310 rows remaining


# Check numbers after filters
table(lab.data$confirmado, useNA = "always")
# Confirmado   Negativo       <NA> 
#   7984623    9693687          0 


# Define Outcome: Create vars about hospitalization and outcome (death)
lab.data <- lab.data %>% mutate.(diff_sample_hosp=dt_interna_dt-dt_coleta,
                        diff_sample_death=dt_evoluca-dt_coleta)
# Define related hospitalization or death: Tests after or before 28 days of hospitalization and death up to 90 days after test.
lab.data <- lab.data %>% mutate.(hosp_event=case_when.(abs(diff_sample_hosp)<=28~"Hospitalization",
                                             TRUE~"No hospitalization"),
                        death_event=case_when.(evolucao %in% c("2","3") &(diff_sample_death>=0 &diff_sample_death<=28)~"Death-related",
                                              TRUE~"No related"),
                        death_event_sensi=case_when.(evolucao %in% c("2","3") &(diff_sample_death>=0 &diff_sample_death<=90)~"Death-related",
                                               TRUE~"No related"))



lab.data <- lab.data %>% drop_na(cod_ibge,sexo,idade,dt_coleta) %>% 
  mutate(uf=factor(str_extract(cod_ibge,"^.{2}")),
         outcome=case_when(
           (hosp_event=="Hospitalization" | death_event=="Death-related")~"Hosp_Death",
           TRUE~"No event"),
         raca=fct_recode(factor(raca),
                         White="1",
                         Black="2",
                         Asian="3",
                         Mixed="4",
                         Indigenous="5"
         ))

#drop_na: removed 2,200 rows (<1%), 17,676,110 rows remaining
# mutate: converted 'raca' from integer to factor (0 new NA)
# new variable 'uf' (factor) with 28 unique values and 0% NA
# new variable 'outcome' (character) with 2 unique values and 0% NA
saveRDS(lab.data,"lab.data_inter0812.rds")


lab.data <- lab.data %>% 
  arrange.(id_vigvac,dt_coleta) %>%
  mutate.(DaysDiffTest = dt_coleta - lags.(dt_coleta),.by = id_vigvac) 

lab.data <- lab.data %>% group_by(id_vigvac)%>% filter(case_when(n()>1 & outcome=="Hosp_Death" & lag(confirmado)=="Confirmado" & confirmado=="Negativo"~DaysDiffTest>28 | is.na(DaysDiffTest),
                   n()>1& outcome=="Hosp_Death"  & lead(confirmado)=="Confirmado" & confirmado=="Negativo"~lead(DaysDiffTest)>28,
                   TRUE~TRUE)) 
# group_by: one grouping variable (id_vigvac)
# filter (grouped): removed 19,243 rows (<1%), 17,656,867 rows remaining
saveRDS(lab.data,"lab.data0812.rds")
toc()

vacc <- arrow::read_parquet("/dados/longs/vacc_clean.parquet")
vacc <- vacc %>% mutate(flag_incon=case_when(
  d2_data_aplicacao-d1_data_aplicacao<13~1L,
  d1_data_aplicacao<"2021-01-18"~1L,
  !is.na(d2_nome_vacina )& d2_nome_vacina!=d1_nome_vacina~1L,
  TRUE~0L
))
vac_test <- left_join.(lab.data,vacc,by="id_vigvac")
vac_test <- vac_test %>% mutate(flag_incon=if_else(is.na(flag_incon),0L,flag_incon))
#mutate: changed 1,249,226 values (7%) of 'flag_incon' (1249226 fewer NA)
vac_test_clean <- vac_test %>% filter(flag_incon==0)
#filter: removed 395,947 rows (2%), 17,260,920 rows remaining
# Recode vaccination
vac_test_clean <- vac_test_clean %>% mutate.(d1_nome_vacina=case_when(
  d1_nome_vacina=="Covid-19-AstraZeneca"~"AZ",
  d1_nome_vacina=="Covid-19-Coronavac-Sinovac/Butantan"~"CV",
  d1_nome_vacina=="Vacina covid-19 - Ad26.COV2.S - Janssen-Cilag"~"Ad26",
  d1_nome_vacina=="Vacina covid-19 - BNT162b2 - BioNTech/Fosun Pharma/Pfizer"~"BNT162b2",
  is.na(d1_nome_vacina)~"uv"),
  d2_nome_vacina=case_when(
    d2_nome_vacina=="Covid-19-AstraZeneca"~"AZ",
    d2_nome_vacina=="Covid-19-Coronavac-Sinovac/Butantan"~"CV",
    d2_nome_vacina=="Vacina covid-19 - Ad26.COV2.S - Janssen-Cilag"~"Ad26",
    d2_nome_vacina=="Vacina covid-19 - BNT162b2 - BioNTech/Fosun Pharma/Pfizer"~"BNT162b2"),
  d3_nome_vacina=case_when(
    d3_nome_vacina=="Covid-19-AstraZeneca"~"AZ",
    d3_nome_vacina=="Covid-19-Coronavac-Sinovac/Butantan"~"CV",
    d3_nome_vacina=="Vacina covid-19 - Ad26.COV2.S - Janssen-Cilag"~"Ad26",
    d3_nome_vacina=="Vacina covid-19 - BNT162b2 - BioNTech/Fosun Pharma/Pfizer"~"BNT162b2")) 


## Change class of variables/create number of comorb
vac_test_clean <- vac_test_clean %>%
  mutate(
    sexo = as.factor(sexo),
    confirmado = as.factor(confirmado),
    cod_ibge = as.factor(cod_ibge),
    tipo_teste = as.factor(tipo_teste),
    days_first_positive=dt_coleta - firstcasedate
  ) %>%
  mutate_rowwise.(n_comorb = sum(c_across.(diabetes:drc))) %>%
  mutate(across(gestante:drc, as.factor)) %>%
  mutate(prev_infected = case_when(
    days_first_positive > 0 ~ "Prev-Infected",
    TRUE ~ "No-Prev-Infected"
  )) %>%
  mutate(
    n_comorb = fct_lump(factor(n_comorb), n = 2, other_level = "2+"),
    days_vacc_from_firstpos = as.numeric(d1_data_aplicacao - firstcasedate),
    raca = fct_explicit_na(factor(raca))
  )






#create days since test
vac_test_clean <- vac_test_clean %>% mutate(days_vacc_1_test = as.numeric(dt_inicio_sintomas - d1_data_aplicacao),
                                            days_vacc_2_test = as.numeric(dt_inicio_sintomas - d2_data_aplicacao),
                                            days_vacc_3_test = as.numeric(dt_inicio_sintomas - d3_data_aplicacao),
                                            vacc_status = case_when(!is.na(d3_data_aplicacao) & days_vacc_3_test>=0 & days_vacc_3_test <=6~"v3_0:6",
                                                                    !is.na(d3_data_aplicacao) & days_vacc_3_test>=7 & days_vacc_3_test <=13~"v3_7:13",
                                                                    !is.na(d3_data_aplicacao) & days_vacc_3_test>=14 & days_vacc_3_test <=30~"v3_14:30",
                                                                    !is.na(d3_data_aplicacao) & days_vacc_3_test>=31~"v3_31+",
                                                                    is.na(d2_data_aplicacao) &  days_vacc_1_test <0 ~ "uv",
                                                                    is.na(d2_data_aplicacao) & days_vacc_1_test >=0 & days_vacc_1_test <= 6 ~ "v1_0:6",
                                                                    is.na(d2_data_aplicacao) & days_vacc_1_test >=7 & days_vacc_1_test <= 13 ~ "v1_7:13",
                                                                    is.na(d2_data_aplicacao) & days_vacc_1_test >=14 ~ "v1_14+",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_1_test <0 ~ "test_before_vaccination1",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test <0 ~ "test_before_vaccination2",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=0 & days_vacc_2_test <=13 ~ "v2_0:13",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=14 & days_vacc_2_test <=30 ~ "v2_14:30",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=31 & days_vacc_2_test <=60 ~ "v2_31:60",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=61 & days_vacc_2_test <=90 ~ "v2_61:90",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=91 & days_vacc_2_test <=120 ~ "v2_91:120",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=121 & days_vacc_2_test <=150 ~ "v2_121:150",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=151 & days_vacc_2_test <=180 ~ "v2_151:180",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >180 ~ "v2_181+"))

#### People vaccinated after test- as unvaccinated
vac_test_clean <- vac_test_clean %>% mutate(vacc_status = case_when(vacc_status == "test_before_vaccination1" ~ "uv",
                                                                  vacc_status == "test_before_vaccination2" & days_vacc_1_test >= 0 & days_vacc_1_test <=6 ~ "v1_0:6",
                                                                  vacc_status == "test_before_vaccination2" & days_vacc_1_test >= 7 & days_vacc_1_test <=13 ~ "v1_7:13",
                                                                  vacc_status == "test_before_vaccination2" & days_vacc_1_test >=14 ~ "v1_14+",
                                                                  TRUE ~ vacc_status))
#Add vaccine type 
vac_test_clean <- vac_test_clean %>% mutate(vs_type = paste(d1_nome_vacina, vacc_status, sep="_"))

#Assign the AZ_uv, CV_uv, BNT162b2_uv and uv_uv to unvaccinated \ Remove Janssen with second dose (inconsistent) 
vac_test_clean <- vac_test_clean %>% mutate(vs_type = case_when(vs_type == "Ad26_uv" ~ "uv",
                                                              vs_type == "CV_uv" ~ "uv",
                                                              vs_type == "BNT162b2_uv" ~ "uv",
                                                              vs_type == "AZ_uv" ~ "uv",
                                                              vs_type == "uv_NA" ~ "uv",
                                                              TRUE ~ vs_type),
                                          vs_type=fct_relevel(factor(vs_type),"uv")) %>% filter(!str_detect(as.character(vs_type),"Ad26_v2")) %>% droplevels()

# mutate: converted 'vs_type' from character to factor (0 new NA)
#filter: removed 5 rows (<1%), 17,260,915 rows remaining

vac_test_clean <- vac_test_clean %>% mutate(vacc_3_type=case_when(
  str_detect(vacc_status,"v3_")~ paste(d1_nome_vacina,d3_nome_vacina,sep = "-")))





# Drop inconsistencies -filter: 
vac_test_clean <- vac_test_clean %>% filter(dt_coleta<"2021-12-09", dt_coleta+21>dt_inicio_sintomas)
#filter: removed 85,189 rows (<1%), 17,175,726 rows remaining
##Create variable for propensity scores-filter
vac_test_clean <- vac_test_clean %>% mutate(psm_vacc=if_else(vacc_status %in% c("uv","v1_0:6","v1_7:13"),"UVacc","Vacc"),
                                            date_year=as.numeric(dt_inicio_sintomas-as.Date("2021-01-01")),
                                            outcome_confirmado=case_when(outcome=="Hosp_Death" & confirmado=="Confirmado"~"Hosp_Death_Positive",
                                                                         outcome!="Hosp_Death" & confirmado=="Confirmado"~"Outpatient_Positive",
                                                                         confirmado=="Negativo"~"Negative")) %>% filter(sexo!="I",sexo!="NA",!is.na(sexo),uf!="99") %>% droplevels()
#filter: removed 1,153 rows (<1%), 17,174,573 rows remaining

# Drop missing values-  ## Recode gestante and puerpera (clear errors), recode age greater or equal 100 as 100
vac_test_clean <- vac_test_clean  %>% drop_na(cod_ibge,sexo,idade,dt_coleta) %>% mutate(uf=fct_relevel(uf,"35"),
                                                                                        puerpera=as.factor(if_else(idade>55,"0",
                                                                                                                   as.character(puerpera))),
                                                                                        gestante=as.factor(if_else(idade>55,"0",
                                                                                                                   as.character(gestante))),
                                                                                        vs_type2=case_when(
                                                                                          str_detect(vs_type,"v3_")~paste(vs_type,d3_nome_vacina, sep="_"),
                                                                                          TRUE~as.character(vs_type)
                                                                                        ))%>% mutate(idade=if_else(idade>99,100,idade),
                                                                                                     vs_type2=fct_relevel(factor(vs_type2),"uv"),
                                                                                                     confirmado=fct_relevel(confirmado,"Negativo"))%>% 
  mutate(length=if_else(str_detect(vs_type2,"CV_v3_"),as.Date("2021-12-08")-d3_data_aplicacao,NULL),
         time_2nd_3rd=if_else(str_detect(vs_type2,"CV_v3_"),d3_data_aplicacao-d2_data_aplicacao,NULL))


#Include microregion code
regiao_imediata <- readxl::read_excel("regioes_geograficas_composicao_por_municipios_2017_20180911.xlsx")
#convert to 6 digit ibge code
regiao_imediata <- regiao_imediata %>% mutate(CD_GEOCODI=str_extract(CD_GEOCODI,"[:digit:]{6}"))
regiao_imediata <- regiao_imediata %>% select(CD_GEOCODI,cod_rgi)
vac_test_clean <- vac_test_clean %>% left_join(regiao_imediata,by=c("cod_ibge"="CD_GEOCODI"))
# Check missing
vac_test_clean %>% filter(is.na(cod_rgi)) %>% count(uf)
#All missing are from DF- assign value of Df 530001 and convert to factor
vac_test_clean <- vac_test_clean %>% mutate(cod_rgi=factor(if_else(is.na(cod_rgi),"530001",cod_rgi)))
vac_test_clean %>% filter(is.na(cod_rgi)) %>% count(uf)

saveRDS(vac_test_clean,"db_tnd_analysis.RDS")



### addendum 
vac_test_clean <- read_rds("db_tnd_analysis.RDS")



#Create new categories
vac_test_clean <- vac_test_clean %>% mutate(days_vacc_1_test = as.numeric(dt_inicio_sintomas - d1_data_aplicacao),
                                            days_vacc_2_test = as.numeric(dt_inicio_sintomas - d2_data_aplicacao),
                                            days_vacc_3_test = as.numeric(dt_inicio_sintomas - d3_data_aplicacao),
                                            vacc_status = case_when(!is.na(d3_data_aplicacao) & days_vacc_3_test>=0 & days_vacc_3_test <=6~"v3_0:6",
                                                                    !is.na(d3_data_aplicacao) & days_vacc_3_test>=7 & days_vacc_3_test <=13~"v3_7:13",
                                                                    !is.na(d3_data_aplicacao) & days_vacc_3_test>=14 & days_vacc_3_test <=30~"v3_14:30",
                                                                    !is.na(d3_data_aplicacao) & days_vacc_3_test>=31~"v3_31+",
                                                                    is.na(d2_data_aplicacao) &  days_vacc_1_test <0 ~ "uv",
                                                                    is.na(d2_data_aplicacao) & days_vacc_1_test >=0 & days_vacc_1_test <= 6 ~ "v1_0:6",
                                                                    is.na(d2_data_aplicacao) & days_vacc_1_test >=7 & days_vacc_1_test <= 13 ~ "v1_7:13",
                                                                    is.na(d2_data_aplicacao) & days_vacc_1_test >=14 ~ "v1_14+",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_1_test <0 ~ "test_before_vaccination1",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test <0 ~ "test_before_vaccination2",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=0 & days_vacc_2_test <=13 ~ "v2_0:13",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=14 & days_vacc_2_test <=30 ~ "v2_14:30",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=31 & days_vacc_2_test <=60 ~ "v2_31:60",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=61 & days_vacc_2_test <=90 ~ "v2_61:90",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=91 & days_vacc_2_test <=120 ~ "v2_91:120",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=121 & days_vacc_2_test <=150 ~ "v2_121:150",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >=151 & days_vacc_2_test <=180 ~ "v2_151:180",
                                                                    !is.na(d2_data_aplicacao) & days_vacc_2_test >180 ~ "v2_181+"))

#### People vaccinated after test- as unvaccinated
vac_test_clean <- vac_test_clean %>% mutate(vacc_status = case_when(vacc_status == "test_before_vaccination1" ~ "uv",
                                                                    vacc_status == "test_before_vaccination2" & days_vacc_1_test >= 0 & days_vacc_1_test <=6 ~ "v1_0:6",
                                                                    vacc_status == "test_before_vaccination2" & days_vacc_1_test >= 7 & days_vacc_1_test <=13 ~ "v1_7:13",
                                                                    vacc_status == "test_before_vaccination2" & days_vacc_1_test >=14 ~ "v1_14+",
                                                                    TRUE ~ vacc_status))
#Add vaccine type 
vac_test_clean <- vac_test_clean %>% mutate(vs_type = paste(d1_nome_vacina, vacc_status, sep="_"))

#Assign the AZ_uv, CV_uv, BNT162b2_uv and uv_uv to unvaccinated \ Remove Janssen with second dose (inconsistent) 
vac_test_clean <- vac_test_clean %>% mutate(vs_type = case_when(vs_type == "Ad26_uv" ~ "uv",
                                                                vs_type == "CV_uv" ~ "uv",
                                                                vs_type == "BNT162b2_uv" ~ "uv",
                                                                vs_type == "AZ_uv" ~ "uv",
                                                                vs_type == "uv_NA" ~ "uv",
                                                                TRUE ~ vs_type),
                                            vs_type=fct_relevel(factor(vs_type),"uv")) %>% filter(!str_detect(as.character(vs_type),"Ad26_v2")) %>% droplevels()

vac_test_clean <- vac_test_clean %>% mutate(vacc_3_type=case_when(
  str_detect(vacc_status,"v3_")~ paste(d1_nome_vacina,d3_nome_vacina,sep = "-")))

vac_test_clean <- vac_test_clean  %>% drop_na(cod_ibge,sexo,idade,dt_coleta) %>% mutate(vs_type2=case_when(
                                                                                          str_detect(vs_type,"v3_")~paste(vs_type,d3_nome_vacina, sep="_"),
                                                                                          TRUE~as.character(vs_type)
                                                                                        ))%>% mutate(length=if_else(str_detect(vs_type2,"CV_v3_"),as.Date("2021-11-11")-d3_data_aplicacao,NULL),
         time_2nd_3rd=if_else(str_detect(vs_type2,"CV_v3_"),d3_data_aplicacao-d2_data_aplicacao,NULL),
         vs_type2=fct_relevel(vs_type2,"uv"))

saveRDS(vac_test_clean,"db_tnd_analysis_rev.RDS")

View = function(x) { utils::View(x)} #make possible to cancel view command
