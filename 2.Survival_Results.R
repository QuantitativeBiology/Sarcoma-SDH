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

clinical_data_clean$SDHB <- as.numeric(RNA$SDHB)
clinical_data_clean$SDHC <- as.numeric(RNA$SDHC)
clinical_data_clean$SDHD <- as.numeric(RNA$SDHD)

clinical_data_clean <- clinical_data_clean[clinical_data_clean$`Sarcoma Histopathological Subtype` == "UPS",]

clinical_data_clean$SDHB <- ifelse(clinical_data_clean$SDHB > mean(clinical_data_clean$SDHB), "High", "Low")

km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~ SDHB, data=clinical_data_clean)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival SDHD")

x



