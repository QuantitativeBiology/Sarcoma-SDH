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


km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~ SDHB, data=clinical_data_clean)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival SDHB")
x

clinical_data_clean$`Local Recurrence` <- as.character(clinical_data_clean$`Local Recurrence`)

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~  Gender + `Sarcoma Histopathological Subtype` + SDHB +  `Neo Adjuvant / Adjuvant Treatment` + `Local Recurrence`  , data = clinical_data_clean)
ggforest(cox)


clinical_data_clean_ups <- clinical_data_clean[clinical_data_clean$`Sarcoma Histopathological Subtype` == "UPS",]

clinical_data_clean_ups$SDHB <- ifelse(clinical_data_clean_ups$SDHB_val > mean(clinical_data_clean_ups$SDHB_val), "High", "Low")

km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~ SDHB, data=clinical_data_clean_ups)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Free Survival SDHB: UPS")
x



