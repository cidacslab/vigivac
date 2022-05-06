##############################################################################
# Name of file: 02_Modelling TND
# Original author(s): Thiago Cerqueira
# Original date: 11 Nov 21
# Latest update author (if not using version control) - thiago.c.silva@fiocruz.br
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 4.1.1
# Description of content: Modelling
# Approximate run time: Unknown
##############################################################################
setwd("/dados/Analysis/Thiago/Scripts/CoronaVac-Waning/TND")
pacman::p_load(tidyverse,tidylog,tidytable,mgcv,gtsummary,rms, tictoc)


# Define functions
models_fit <- function(dataset){
  if(min(dataset$idade)<56)
    bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
        family = binomial, nthreads = 24, discrete = T,data = dataset)
  else
    bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb,
        family = binomial, nthreads = 24, discrete = T,data = dataset)
}

models_check <- function(dataset){
    mod_minimal <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ sexo,
        family = binomial, nthreads = 24, discrete = T,data = dataset)
    values_1 <- as.data.frame(summary(mod_minimal)$p.table) %>%
      rownames_to_column(var = "term") %>%
      select(1:3) %>%
      filter(str_detect(term, "vs_type2CV")) %>% 
      filter(case_when(str_detect(term,"_v3_")~str_detect(term,"BNT"),
                       TRUE~TRUE)) %>% 
      mutate(mod = "minimal")
    mod_mod <-bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo,
                       family = binomial, nthreads = 24, discrete = T,data = dataset)
    values_2 <- as.data.frame(summary(mod_mod)$p.table) %>%
      rownames_to_column(var = "term") %>%
      select(1:3) %>%
      filter(str_detect(term, "vs_type2CV")) %>%
      filter(case_when(str_detect(term,"_v3_")~str_detect(term,"BNT"),
                       TRUE~TRUE)) %>% 
      mutate(mod = "medium")
    mod_ms <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante,
                      family = binomial, nthreads = 24, discrete = T,data = dataset)
    values_3 <- as.data.frame(summary(mod_ms)$p.table) %>%
      rownames_to_column(var = "term") %>%
      select(1:3) %>%
      filter(str_detect(term, "vs_type2CV")) %>%
      filter(case_when(str_detect(term,"_v3_")~str_detect(term,"BNT"),
                       TRUE~TRUE)) %>% 
      mutate(mod = "ms")
    mod_full <- bam(confirmado=="Confirmado" ~ s(date_year,bs="cr")+s(idade,bs="cr")+vs_type2+ prev_infected+ sexo +uf +n_comorb+ puerpera + gestante+raca,
                  family = binomial, nthreads = 24, discrete = T,data = dataset)
    values_4 <- as.data.frame(summary(mod_full)$p.table) %>%
      rownames_to_column(var = "term") %>%
      select(1:3) %>%
      filter(str_detect(term, "vs_type2CV")) %>%
      filter(case_when(str_detect(term,"_v3_")~str_detect(term,"BNT"),
                       TRUE~TRUE)) %>% 
      mutate(mod = "full")
    results <- bind_rows(values_1,values_2,values_3,values_4)
    return(results)
}

extract_ve <- function(model, name, booster) {
  round(exp(cbind("Odds Ratio" = coef(model), confint.default(model, level = 0.95))), digits = 3) %>%
    as.data.frame() %>%
    rownames_to_column(var = "term") %>%
    filter(str_detect(term, name)) %>%
    filter(case_when(str_detect(term,"_v3_")~str_detect(term,booster),
                     TRUE~TRUE)) %>% 
    mutate(term = str_remove(term, paste0(name,"_"))) %>%
    mutate(across(c(2:4), ~ (1 - .x) * 100, .names = "{.col}_ve")) %>%
    select(1, 5:7) %>%
    mutate(ve_ci = paste0(sprintf("%.1f", `Odds Ratio_ve`), "% (",sprintf("%.1f",`97.5 %_ve`), "-",sprintf("%.1f",`2.5 %_ve`),")")) %>% 
    mutate(term=factor(term,levels=c("v1_0:6","v1_7:13","v1_14+","v2_0:13",
                                     "v2_14:30","v2_31:60","v2_61:90","v2_91:120",
                                     "v2_121:150","v2_151:180","v2_181+","v3_0:6_BNT162b2",
                                     "v3_7:13_BNT162b2","v3_14:30_BNT162b2","v3_31+_BNT162b2"))) %>% arrange(term) 
}

select(confirmado, vs_type2) %>% filter(str_detect(vs_type2, "CV_|uv")) %>% 
  filter(case_when(
    str_detect(vs_type2, "_v3_") ~ str_detect(vs_type2, "BNT"),
    TRUE ~ TRUE
  )) %>% droplevels() %>% tbl_summary(by = confirmado, percent = "row") %>% add_overall(last = T) 

extract_numbers_overall <- function(data){
  data %>%
    select(confirmado, vs_type2) %>%
    mutate(vs_type2 = fct_relevel(
      vs_type2, "uv", "CV_v1_0:6",
      "CV_v1_7:13", "CV_v1_14+", "CV_v2_0:13",
      "CV_v2_14:30", "CV_v2_31:60", "CV_v2_61:90",
      "CV_v2_91:120", "CV_v2_121:150", "CV_v2_151:180",
      "CV_v2_181+", "CV_v3_0:6_BNT162b2", "CV_v3_7:13_BNT162b2", "CV_v3_14:30_BNT162b2", "CV_v3_31+_BNT162b2"
    )) %>%
    filter(case_when(
      str_detect(vs_type2, "_v3_") ~ str_detect(vs_type2, "BNT"),
      TRUE ~ TRUE
    )) %>%
    filter(str_detect(vs_type2, "CV_|uv")) %>%
    droplevels() %>%
    tbl_summary(by = confirmado, percent = "row", digits = list(all_categorical() ~ c(0, 1))) %>%
    add_overall(last = T) 
}

extract_numbers_cat <- function(data){
  data %>%
    nest(-idade_cat) %>%
    mutate(events = map(data, . %>% select(confirmado, vs_type2) %>% mutate(vs_type2 = fct_relevel(
      vs_type2, "uv", "CV_v1_0:6",
      "CV_v1_7:13", "CV_v1_14+", "CV_v2_0:13",
      "CV_v2_14:30", "CV_v2_31:60", "CV_v2_61:90",
      "CV_v2_91:120", "CV_v2_121:150", "CV_v2_151:180",
      "CV_v2_181+", "CV_v3_0:6_BNT162b2", "CV_v3_7:13_BNT162b2", "CV_v3_14:30_BNT162b2", "CV_v3_31+_BNT162b2"
    )) %>% filter(case_when(
      str_detect(vs_type2, "_v3_") ~ str_detect(vs_type2, "BNT"),
      TRUE ~ TRUE
    )) %>% filter(str_detect(vs_type2, "CV_|uv")) %>% 
      droplevels() %>% tbl_summary(by = confirmado, percent = "row", digits = list(all_categorical() ~ c(0, 1))) %>% 
      add_overall(last = T))) %>%
    select(events)
}
##############################################################################

vac_test_clean <- read_rds("db_tnd_analysis_rev.RDS")

# Datasets 
# reorder/create cat
vac_test <- vac_test_clean %>% mutate(idade_10=cut(idade, breaks=c(18,seq(30,80,by=10),Inf), right = FALSE,
                                                                                           labels=c(paste(c(18,seq(30,70, by=10)),seq(29,79, by=10),sep="-"),"80+"), ordered_result=T,
                                                                                           include.lowest=T),
                                                                      idade_cat=cut(idade, breaks=c(18,60,80,+Inf), right = FALSE,
                                                                                   labels=c(paste(c(18,60),c(59,79),sep="-"),"80+"), ordered_result=T,
                                                                                   include.lowest=T),
                                      vs_type2=fct_relevel(vs_type2,"uv"))


## Infection
vac_test_rt_pcr <- vac_test %>% filter(tipo_teste=="RT-PCR")

## Severe
vac_test_rt_pcr_severe <- vac_test%>% filter(tipo_teste=="RT-PCR") %>% filter(outcome_confirmado!="Outpatient_Positive") %>% droplevels()

vac_seve_antigen <- vac_test%>% filter(outcome_confirmado!="Outpatient_Positive") %>% droplevels()


## VE- Models
overall_db <- tibble(outcome=c("Infection","Severe"), data=list(vac_test_rt_pcr,vac_test_rt_pcr_severe))
overall_ve <- overall_db  %>%  mutate(fit = map(data,models_fit ))
overall_ve <- overall_ve%>% mutate(ve=map(fit,extract_ve,name="vs_type2CV",booster="BNT"))
overall_ve %>% select(outcome,ve) %>% unnest(ve) %>% View()

## by age-  infection
vac_df_cat_infection <- vac_test_rt_pcr %>% nest(-idade_cat) %>%  mutate(fit = map(data,models_fit ))
vac_df_cat_infection %>% mutate(ve=map(fit,extract_ve,name="vs_type2CV",booster="BNT"))  %>%  unnest(ve) %>% select(-c(data,fit)) %>% View()


## by age-  severe
vac_df_cat_severe <- vac_test_rt_pcr_severe %>% nest(-idade_cat) %>%  mutate(fit = map(data,models_fit ))
vac_df_cat_severe%>% mutate(ve=map(fit,extract_ve,name="vs_type2CV",booster="BNT"))  %>%  unnest(ve) %>% select(-c(data,fit)) %>% View()


# Check confounders
check_reviewer <- models_check(vac_test_rt_pcr)

# Individuals outcomes

## Hosp

vac_hosp <- vac_test_rt_pcr_severe %>%
  filter(hosp_event == "Hospitalization" & confirmado == "Confirmado" | confirmado == "Negativo") %>%
  nest(-idade_cat) %>%
  mutate(fit = map(data, models_fit))
vac_hosp %>%
  mutate(ve = map(fit, extract_ve, name = "vs_type2CV", booster = "BNT")) %>%
  unnest(ve) %>%
  select(-c(data, fit)) %>%
  View()

## Death

vac_death <- vac_test_rt_pcr_severe %>%
  filter(death_event == "Death-related" & confirmado == "Confirmado" | confirmado == "Negativo") %>%
  nest(-idade_cat) %>%
  mutate(fit = map(data, models_fit))
vac_death %>%
  mutate(ve = map(fit, extract_ve, name = "vs_type2CV", booster = "BNT")) %>%
  unnest(ve) %>%
  select(-c(data, fit)) %>%
  View()


death_db <- vac_test_rt_pcr_severe %>% filter(death_event == "Death-related" & confirmado == "Confirmado" | confirmado == "Negativo")
hosp_db <- vac_test_rt_pcr_severe %>% filter(hosp_event == "Hospitalization" & confirmado == "Confirmado" | confirmado == "Negativo")
hospdeath_db <- tibble(outcome = c("Death", "Hosp"), data = list(death_db, hosp_db))
hospdeath_ve <- hospdeath_db %>% mutate(fit = map(data, models_fit))
hospdeath_ve <- hospdeath_ve %>%
  mutate(ve = map(fit, extract_ve, name = "vs_type2CV", booster = "BNT")) %>%
  unnest(ve) %>%
  select(-c(data, fit))



#Numbers
#Infection
vac_test_rt_pcr %>% extract_numbers_overall()
no_rt <- vac_test_rt_pcr %>% extract_numbers_cat()
ag_pcr <- vac_test %>% extract_numbers_overall()
#Overall-severe
vac_test_rt_pcr_severe %>% extract_numbers_overall()

vac_seve_antigen%>% extract_numbers_overall()
#severe
number_severe <- vac_test_rt_pcr_severe %>% extract_numbers_cat()


#only hosp
number_hosp <- vac_test_rt_pcr_severe %>%
  filter(hosp_event=="Hospitalization" & confirmado=="Confirmado" | confirmado=="Negativo") %>% 
  nest(-idade_cat) %>%
  mutate(events = map(data, . %>% select(confirmado, vs_type2) %>% filter(case_when(
    str_detect(vs_type2, "_v3_") ~ str_detect(vs_type2, "BNT"),
    TRUE ~ TRUE
  )) %>% filter(str_detect(vs_type2, "CV_|uv")) %>% droplevels() %>% tbl_summary(by = confirmado) %>% add_overall(last = T))) %>%
  select(events)

#only death
number_death <- vac_test_rt_pcr_severe %>%
  filter(death_event=="Death-related" & confirmado=="Confirmado" | confirmado=="Negativo") %>% 
  nest(-idade_cat) %>%
  mutate(events = map(data, . %>% select(confirmado, vs_type2) %>% filter(case_when(
    str_detect(vs_type2, "_v3_") ~ str_detect(vs_type2, "BNT"),
    TRUE ~ TRUE
  )) %>% filter(str_detect(vs_type2, "CV_|uv")) %>% droplevels() %>% tbl_summary(by = confirmado) %>% add_overall(last = T))) %>%
  select(events)


# Change reference category
vac_test_rt_pcr_seni <- vac_test %>%
  filter(tipo_teste == "RT-PCR") %>%
  mutate(vs_type2 = fct_relevel(vs_type2, "CV_v2_181+"))

vac_rt_pcr_severe_sensi <- vac_test %>%
  filter(tipo_teste == "RT-PCR") %>%
  filter(outcome_confirmado != "Outpatient_Positive") %>%
  mutate(vs_type2 = fct_relevel(vs_type2, "CV_v2_181+")) %>%
  droplevels()

infection_sensi <- bam(confirmado == "Confirmado" ~ s(date_year, bs = "cr") + s(idade, bs = "cr") + vs_type2 + prev_infected + sexo + uf + n_comorb + puerpera + gestante,
  family = binomial, nthreads = 24, discrete = T, data = vac_test_rt_pcr_seni
)
severe_sensi <- bam(confirmado == "Confirmado" ~ s(date_year, bs = "cr") + s(idade, bs = "cr") + vs_type2 + prev_infected + sexo + uf + n_comorb + puerpera + gestante,
  family = binomial, nthreads = 24, discrete = T, data = vac_rt_pcr_severe_sensi
)

infection_sensi_180_ag <- bam(confirmado == "Confirmado" ~ s(date_year, bs = "cr") + s(idade, bs = "cr") + vs_type2 + prev_infected + sexo + uf + n_comorb + puerpera + gestante,
  family = binomial, nthreads = 24, discrete = T, data = vac_test %>% mutate(vs_type2 = fct_relevel(vs_type2, "CV_v2_181+"))
)

severe_sensi_180_ag <- bam(confirmado == "Confirmado" ~ s(date_year, bs = "cr") + s(idade, bs = "cr") + vs_type2 + prev_infected + sexo + uf + n_comorb + puerpera + gestante,
  family = binomial, nthreads = 24, discrete = T, data = vac_test %>% filter(outcome_confirmado != "Outpatient_Positive") %>% mutate(vs_type2 = fct_relevel(vs_type2, "CV_v2_181+")) %>% droplevels()
)

extract_ve(severe_sensi_180_ag, name = "vs_type2CV", booster = "BNT") %>% View()

infection_cat_sensi <- vac_test_rt_pcr_seni %>%
  nest(-idade_cat) %>%
  mutate(fit = map(data, models_fit)) %>%
  mutate(ve = map(fit, extract_ve, name = "vs_type2CV", booster = "BNT")) %>%
  unnest(ve) %>%
  select(-c(data, fit)) %>%
  filter(str_detect(term, "v3"))

infection_cat_sensi_ag <- vac_test %>% mutate(vs_type2=fct_relevel(vs_type2,"CV_v2_181+")) %>% 
  nest(-idade_cat) %>%
  mutate(fit = map(data, models_fit)) %>%
  mutate(ve = map(fit, extract_ve, name = "vs_type2CV", booster = "BNT")) %>%
  unnest(ve) %>%
  select(-c(data, fit)) %>%
  filter(str_detect(term, "v3"))


severe_cat_sensi <- vac_rt_pcr_severe_sensi %>%
  nest(-idade_cat) %>%
  mutate(fit = map(data, models_fit)) %>%
  mutate(ve = map(fit, extract_ve, name = "vs_type2CV", booster = "BNT")) %>%
  unnest(ve) %>%
  select(-c(data, fit)) %>%
  filter(str_detect(term, "v3"))

severe_cat_sensi_ag <- vac_test%>% filter(outcome_confirmado!="Outpatient_Positive")%>% mutate(vs_type2=fct_relevel(vs_type2,"CV_v2_181+")) %>% droplevels() %>%
  nest(-idade_cat) %>%
  mutate(fit = map(data, models_fit)) %>%
  mutate(ve = map(fit, extract_ve, name = "vs_type2CV", booster = "BNT")) %>%
  unnest(ve) %>%
  select(-c(data, fit)) %>%
  filter(str_detect(term, "v3"))



# Sensitivity -antigen included

antigen_db <- tibble(outcome=c("Infection","Severe"), data=list(vac_test_clean,vac_severe_sensi))
antigen_ve <- antigen_db  %>%  mutate(fit = map(data,models_fit ))
antigen_ve <- antigen_ve%>% mutate(ve=map(fit,extract_ve,name="vs_type2CV",booster="BNT"))
antigen_ve %>% select(outcome,ve) %>% unnest(ve) %>% View()


# Infection

## Tables



vac_test_infection %>%
  group_by(id_vigvac) %>%  mutate(unique_id=row_number(id_vigvac)) %>% ungroup() %>% 
  mutate( vs_type3 = case_when(str_detect(vs_type2, "CV_")~ "CoronaVac",
                         vs_type2=="uv"~"Unvaccinated",
                         TRUE~"Other Vaccines")
  ) %>% select(age_group, vs_type3, confirmado) %>%
  tbl_strata(
    strata = age_group,
    .tbl_fun =
      ~ .x %>%
        tbl_summary(by = confirmado) %>%
        add_n()
  )%>% bold_labels() 


vac_test  %>%  group_by(id_vigvac) %>%  mutate(unique_id=row_number(id_vigvac)) %>% ungroup() %>% 
  group_by(id_vigvac) %>%  mutate(unique_id=row_number(id_vigvac)) %>% ungroup() %>% 
  mutate( vs_type3 = case_when(str_detect(vs_type2, "CV_")~ "CoronaVac",
                               vs_type2=="uv"~"Unvaccinated",
                               TRUE~"Other Vaccines")
  ,vs_type3=case_when(
      vs_type2=="_CV"~"Other Vaccines",
      vs_type2=="_Ad26"~"Other Vaccines",
      vs_type2=="_AZ"~"Other Vaccines",
      TRUE~vs_type3
    )) %>%
  mutate(
    macro_region = str_extract(cod_ibge, "^.{1}"),
    macro_region = case_when(
      macro_region == "1" ~ "North",
      macro_region == "2" ~ "Northeast",
      macro_region == "3" ~ "Southeast",
      macro_region == "4" ~ "Central-west",
      macro_region == "5" ~ "South"
    )
  ) %>%
  select(idade, sexo, raca,unique_id, idade_cat,tipo_teste, macro_region, gestante, puerpera, diabetes, obesidade, imunossupressao, dcardiaca, drc, n_comorb, prev_infected, hosp_event, death_event,outcome,tipo_teste, vs_type3, confirmado) %>%
  tbl_summary(by = confirmado,digits = list(all_categorical() ~ c(0, 1)),value = list(gestante:drc ~ "1", unique_id~"1"))%>% 
  bold_labels() %>% add_overall(last=T) %>% add_n()




vac_test_infection%>%  filter(str_detect(vs_type2,"CV_v3")) %>% droplevels() %>% 
  mutate(
    age_group = case_when(
      idade < 60 ~ "18-59",
      idade < 80 ~ "60-79",
      idade > 79 ~ "80+"
    ),
    vs_type3 = case_when(
      str_detect(vs_type2, "CV_") ~ as.character(vs_type2),
      vs_type2 == "uv" ~ "Unvaccinated",
      TRUE ~ "Other Vaccines"
    ),
    vs_type3=case_when(
      vs_type3=="CV_v3_CV"~"Other Vaccines",
      vs_type3=="CV_v3_Ad26"~"Other Vaccines",
      vs_type3=="CV_v3_AZ"~"Other Vaccines",
      TRUE~vs_type3
    )
  ) %>%
  mutate(
    macro_region = str_extract(cod_ibge, "^.{1}"),
    macro_region = case_when(
      macro_region == "1" ~ "North",
      macro_region == "2" ~ "Northeast",
      macro_region == "3" ~ "Southeast",
      macro_region == "4" ~ "Central-west",
      macro_region == "5" ~ "South"
    )
  )  %>% filter(vs_type2=="CV_v3_BNT162b2") %>% 
  group_by(id_vigvac) %>%  mutate(unique_id=row_number(id_vigvac)) %>% ungroup() %>% 
  select(age_group,vs_type3, length,idade, sexo, raca,unique_id, age_group,tipo_teste, macro_region, gestante, puerpera, diabetes, obesidade, imunossupressao, dcardiaca, drc, n_comorb, prev_infected, hosp_event, death_event,outcome,tipo_teste, vs_type3, confirmado) %>%
      tbl_summary(by = vs_type3,digits = list(all_categorical() ~ c(0, 1)))%>% add_overall(last=T)%>% bold_labels() 




# Severe
vac_test_severe_sensitivity <- vac_test_model %>% mutate(idade=if_else(idade>99,100L,idade))%>% filter(outcome_confirmado!="Outpatient_Positive") %>% droplevels()
vac_test_severe <- vac_test_model %>% filter(outcome_confirmado!="Outpatient_Positive", tipo_teste=="RT-PCR") %>% droplevels()







vac_test_model <- vac_test_model %>% mutate(vaccinated_unvacc=if_else(vs_type=="uv","uv","v"),
                                          cv_other=case_when(
                                            str_detect(vs_type,"CV")~"CoronaVac",
                                            vs_type!="uv"~ "other vaccines"
                                          ))

vac_test_infection <- vac_test_sensitivity %>% mutate(vaccinated_unvacc=if_else(vs_type=="uv","uv","v"),
                                            cv_other=case_when(
                                              str_detect(vs_type,"CV")~"CoronaVac",
                                              vs_type!="uv"~ "other vaccines",
                                              vs_type=="uv"~"Unvaccinated"
                                            ),
                                            cv_time=case_when(
                                              str_detect(vs_type,"CV")~as.character(vs_type),
                                              !str_detect(vs_type,"CV") & vs_type!="uv" ~ "other vaccines",
                                              vs_type=="uv"~"Unvaccinated"
                                            ))

vac_test_clean %>%
  select(idade, sexo, gestante, puerpera, diabetes, obesidade, imunossupressao, dcardiaca, drc, raca, n_comorb, prev_infected, tipo_teste, confirmado,hosp_event,death_event,outcome) %>%
  tbl_summary(by = confirmado) %>%
  bold_labels()


vac_test_infection %>%
  select(cv_time,confirmado) %>%
  tbl_summary(by = confirmado) %>%
  bold_labels() %>% add_overall(last=T)



#Old approach
models_infection <- list(m1_infection,m1_infection_sensitivity,m1_symp, m18_infection,m60_infection,m80_infection)
names(models_infection) <- c("Infection","Sensitivity", "Symptomatic","18-59","60-79","80+")
map_dfr(models_infection,extract_ve,.id="Model") %>% 
  mutate(term=factor(term,levels=c("v1_0:6","v1_7:13","v1_14+","v2_0:13",
                                   "v2_14:30","v2_31:60","v2_61:90","v2_91:120",
                                   "v2_121:150","v2_151:180","v2_181+","v3_0:7_BNT162b2",
                                   "v3_8:13_BNT162b2","v3_14:30_BNT162b2","v3_31+_BNT162b2"))) %>% arrange(Model,term) %>% View()
