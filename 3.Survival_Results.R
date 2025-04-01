library(tidyverse)
library(survival)
library(survminer)
library(readxl)
library(forestploter)


clinical_data_clean <- read.csv("FILES/clinical_data_clean.csv", row.names=1)


colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Sarcoma.Histopathological.Subtype"] <- "Sarcoma Histopathological Subtype"

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Neo.Adjuvant...Adjuvant.Treatment"] <- "Neo Adjuvant / Adjuvant Treatment"

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Local.Recurrence"] <- "Local Recurrence"

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Distant.Recurrence"] <- "Distant Recurrence"

RNA <- read.csv("FILES/normalized_counts.csv", row.names=1)

colnames(RNA) <- sub("^X", "", colnames(RNA))

common_patients <- intersect(colnames(RNA),row.names(clinical_data_clean))

clinical_data_clean <- clinical_data_clean[common_patients,]
RNA <- RNA[,common_patients]

RNA <- data.frame(t(RNA))

clinical_data_clean$SDHA_val <- as.numeric(RNA$SDHA)
clinical_data_clean$SDHB_val <- as.numeric(RNA$SDHB)
clinical_data_clean$SDHC_val <- as.numeric(RNA$SDHC)
clinical_data_clean$SDHD_val <- as.numeric(RNA$SDHD)


#no significant results comparing subtypes
clinical_data_clean$SDHA <- ifelse(clinical_data_clean$SDHA_val > mean(clinical_data_clean$SDHA_val), "High", "Low")
clinical_data_clean$SDHB <- ifelse(clinical_data_clean$SDHB_val > mean(clinical_data_clean$SDHB_val), "High", "Low")
clinical_data_clean$SDHC <- ifelse(clinical_data_clean$SDHC_val > mean(clinical_data_clean$SDHC_val), "High", "Low")
clinical_data_clean$SDHD <- ifelse(clinical_data_clean$SDHD_val > mean(clinical_data_clean$SDHD_val), "High", "Low")


km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~ SDHD, data=clinical_data_clean)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival SDHD")

clinical_data_clean$`Local Recurrence` <- as.character(clinical_data_clean$`Local Recurrence`)


#clinical_data_clean$SDHB <- relevel(clinical_data_clean$SDHB, ref = "Low")

clinical_data_clean$SDHB <- factor(clinical_data_clean$SDHB, levels = c("Low", "High"))

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  Gender + `Sarcoma Histopathological Subtype` + SDHB +  `Neo Adjuvant / Adjuvant Treatment` + `Local Recurrence`  , data = clinical_data_clean)
ggforest(cox, fontsize = 1)

clinical_data_clean_mfs <- clinical_data_clean

clinical_data_clean_mfs$TIME_TO_METASTISIS <- as.Date(clinical_data_clean_mfs$Data.metastização) - as.Date(clinical_data_clean_mfs$Data.cirurgia)

clinical_data_clean_mfs <- clinical_data_clean_mfs[clinical_data_clean_mfs$TIME_TO_METASTISIS > 0,]

clinical_data_clean_mfs$Distant_Recurrence <- clinical_data_clean_mfs$`Distant Recurrence`

km_fit <- survfit(Surv(TIME_TO_METASTISIS, Distant_Recurrence) ~ SDHD, data=clinical_data_clean_mfs)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Metastastis Free Survival SDHD")
x

clinical_data_clean_rfs <- clinical_data_clean

clinical_data_clean_rfs$Local_Recurrence <- clinical_data_clean_rfs$`Local Recurrence`

clinical_data_clean_rfs$Local_Recurrence[is.na(clinical_data_clean_rfs$Local_Recurrence)] <- "0"

clinical_data_clean_rfs$Local_Recurrence <- as.numeric(clinical_data_clean_rfs$Local_Recurrence)

clinical_data_clean_rfs$TIME_TO_LOCAL_RECURRENCE <- as.Date(clinical_data_clean_rfs$Data.recidiva) - as.Date(clinical_data_clean_rfs$Data.cirurgia)

clinical_data_clean_rfs <- clinical_data_clean_rfs[clinical_data_clean_rfs$TIME_TO_LOCAL_RECURRENCE > 0,]

km_fit <- survfit(Surv(TIME_TO_LOCAL_RECURRENCE, Local_Recurrence) ~ SDHD, data=clinical_data_clean_rfs)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Recurrence Free Survival SDHD")
x

# Calculate PFS #

clinical_data_pfs <- clinical_data_clean

clinical_data_pfs$TIME_PFS <- NA
clinical_data_pfs$PFS <- NA

clinical_data_pfs$Local_Recurrence <- clinical_data_pfs$`Local Recurrence`
clinical_data_pfs$Distant_Recurrence <- clinical_data_pfs$`Distant Recurrence`

clinical_data_pfs$TIME_TO_RECIDIVE <- as.numeric(as.Date(clinical_data_pfs$Data.recidiva) - as.Date(clinical_data_pfs$Data.cirurgia))
clinical_data_pfs$TIME_TO_METASTISIS <- as.numeric(as.Date(clinical_data_pfs$Data.metastização) - as.Date(clinical_data_pfs$Data.cirurgia))

clinical_data_pfs$TIME_TO_METASTISIS
clinical_data_pfs <- clinical_data_pfs[clinical_data_pfs$TIME_TO_RECIDIVE >0,]
clinical_data_pfs <- clinical_data_pfs[clinical_data_pfs$TIME_TO_METASTISIS >0,]


for (patient in row.names(clinical_data_pfs)){
  row <- clinical_data_pfs[patient,]
  if (!is.na(row$TIME_TO_RECIDIVE) && !is.na(row$TIME_TO_METASTISIS) && row$TIME_TO_RECIDIVE != row$TIME_TO_METASTISIS) {
    clinical_data_pfs[patient,]$TIME_PFS <- min(row$TIME_TO_RECIDIVE, row$TIME_TO_METASTISIS)
    clinical_data_pfs[patient,]$PFS <- 1
  } else {
    clinical_data_pfs[patient,]$TIME_PFS <- row$TIME_DEATH_FROM_SURGERY
    clinical_data_pfs[patient,]$PFS <- 0
  }
}

km_fit <- survfit(Surv(TIME_PFS, PFS) ~ SDHC, data=clinical_data_pfs)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Progression Free Survival SDHC")
x


## Inside UPS group only

clinical_data_clean_ups <- clinical_data_clean[clinical_data_clean$`Sarcoma Histopathological Subtype` == "UPS",]



clinical_data_clean_ups$SDHA <- ifelse(clinical_data_clean_ups$SDHA_val > mean(clinical_data_clean_ups$SDHA_val), "High", "Low")
clinical_data_clean_ups$SDHB <- ifelse(clinical_data_clean_ups$SDHB_val > mean(clinical_data_clean_ups$SDHB_val), "High", "Low")
clinical_data_clean_ups$SDHC <- ifelse(clinical_data_clean_ups$SDHC_val > mean(clinical_data_clean_ups$SDHC_val), "High", "Low")
clinical_data_clean_ups$SDHD <- ifelse(clinical_data_clean_ups$SDHD_val > mean(clinical_data_clean_ups$SDHD_val), "High", "Low")

km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~ SDHA, data=clinical_data_clean_ups)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Free Survival SDHB: UPS")
x

clinical_data_clean_ups$SDHB <- factor(clinical_data_clean_ups$SDHB, levels = c("Low", "High"))

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  Gender + SDHB +  `Neo Adjuvant / Adjuvant Treatment` + `Local Recurrence`  , data = clinical_data_clean_ups)
ggforest(cox, fontsize = 1)


clinical_data_clean_mfs <- clinical_data_clean_ups

clinical_data_clean_mfs$TIME_TO_METASTISIS <- as.Date(clinical_data_clean_mfs$Data.metastização) - as.Date(clinical_data_clean_mfs$Data.cirurgia)

clinical_data_clean_mfs <- clinical_data_clean_mfs[clinical_data_clean_mfs$TIME_TO_METASTISIS > 0,]

clinical_data_clean_mfs$Distant_Recurrence <- clinical_data_clean_mfs$`Distant Recurrence`

km_fit <- survfit(Surv(TIME_TO_METASTISIS, Distant_Recurrence) ~ SDHB, data=clinical_data_clean_mfs)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Metastastis Free Survival SDHB: UPS")
x


clinical_data_clean_rfs <- clinical_data_clean_ups

clinical_data_clean_rfs$Local_Recurrence <- clinical_data_clean_rfs$`Local Recurrence`

clinical_data_clean_rfs$Local_Recurrence[is.na(clinical_data_clean_rfs$Local_Recurrence)] <- "0"

clinical_data_clean_rfs$Local_Recurrence <- as.numeric(clinical_data_clean_rfs$Local_Recurrence)

clinical_data_clean_rfs$TIME_TO_LOCAL_RECURRENCE <- as.Date(clinical_data_clean_rfs$Data.recidiva) - as.Date(clinical_data_clean_rfs$Data.cirurgia)

clinical_data_clean_rfs <- clinical_data_clean_rfs[clinical_data_clean_rfs$TIME_TO_LOCAL_RECURRENCE > 0,]

km_fit <- survfit(Surv(TIME_TO_LOCAL_RECURRENCE, Local_Recurrence) ~ SDHB, data=clinical_data_clean_rfs)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Recurrence Free Survival SDHB: UPS")
x


clinical_data_pfs <- clinical_data_clean_ups

clinical_data_pfs$TIME_PFS <- NA
clinical_data_pfs$PFS <- NA

clinical_data_pfs$Local_Recurrence <- clinical_data_pfs$`Local Recurrence`
clinical_data_pfs$Distant_Recurrence <- clinical_data_pfs$`Distant Recurrence`

clinical_data_pfs$TIME_TO_RECIDIVE <- as.numeric(as.Date(clinical_data_pfs$Data.recidiva) - as.Date(clinical_data_pfs$Data.cirurgia))
clinical_data_pfs$TIME_TO_METASTISIS <- as.numeric(as.Date(clinical_data_pfs$Data.metastização) - as.Date(clinical_data_pfs$Data.cirurgia))

clinical_data_pfs$TIME_TO_METASTISIS
clinical_data_pfs <- clinical_data_pfs[clinical_data_pfs$TIME_TO_RECIDIVE > 0,]
clinical_data_pfs <- clinical_data_pfs[clinical_data_pfs$TIME_TO_METASTISIS >0,]


for (patient in row.names(clinical_data_pfs)){
  row <- clinical_data_pfs[patient,]
  if (!is.na(row$TIME_TO_RECIDIVE) && !is.na(row$TIME_TO_METASTISIS) && row$TIME_TO_RECIDIVE != row$TIME_TO_METASTISIS) {
    clinical_data_pfs[patient,]$TIME_PFS <- min(row$TIME_TO_RECIDIVE, row$TIME_TO_METASTISIS)
    clinical_data_pfs[patient,]$PFS <- 1
  } else {
    clinical_data_pfs[patient,]$TIME_PFS <- row$TIME_DEATH_FROM_SURGERY
    clinical_data_pfs[patient,]$PFS <- 0
  }
}

km_fit <- survfit(Surv(TIME_PFS, PFS) ~ SDHB, data=clinical_data_pfs)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Progression Free Survival SDHB: UPS")
x

clinical_data_clean_ups <- clinical_data_clean_ups[clinical_data_clean_ups$`Distant Recurrence` == 1,]

clinical_data_clean_ups$TIME_FROM_METASTISIS_DEATH <- as.numeric(as.Date(clinical_data_clean_ups$Data.UF) - as.Date(clinical_data_clean_ups$Data.metastização))

clinical_data_clean_ups$TIME_FROM_METASTISIS_DEATH
clinical_data_clean_ups <- clinical_data_clean_ups[clinical_data_clean_ups$TIME_FROM_METASTISIS_DEATH >0,]

km_fit <- survfit(Surv(TIME_FROM_METASTISIS_DEATH, Morte.S.N) ~ SDHB, data=clinical_data_clean_ups)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival From Metastasis Date SDHB: UPS")
x

clinical_data_clean_ups$SDHB <- factor(clinical_data_clean_ups$SDHB, levels = c("Low", "High"))

cox <- coxph(Surv(TIME_FROM_METASTISIS_DEATH, Morte.S.N) ~  Gender + SDHB +  `Neo Adjuvant / Adjuvant Treatment` + `Local Recurrence`  , data = clinical_data_clean_ups)
ggforest(cox, fontsize = 1)



