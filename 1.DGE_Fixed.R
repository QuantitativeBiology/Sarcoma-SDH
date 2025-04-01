# Load Libraries
library(dplyr)
library(readxl)
library(limma)
library(edgeR)
library(EnhancedVolcano)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggpubr)
library(RColorBrewer)
library(rempsyc)

# Read data files
RNA_data <- read.csv("FILES/22-2290 Roche RNA GEP-COUNTS 01MAY2023.txt", sep = "\t")
DNA_data <- read_excel("FILES/OSPL_n=82_F1CDXRNA-RMC-RET-22-2290_SG44174_21FEB2023145617.xlsx")

# Filter DNA data
DNA_data <- DNA_data[DNA_data$`PASS/QUALIFIED STATUS_DNA` != "Fail", ]

# Preprocess RNA Data
RNA_data <- subset(RNA_data, select = -c(QC_STATUS, BAITSET, TEST_TYPE, EXPRESSION_UNIT, QC_FLAGS))
row.names(RNA_data) <- sub("\\*.*", "", RNA_data$SPECIMEN)
RNA_data$SPECIMEN <- NULL

# Find common FMIs
common_fmis <- intersect(DNA_data$`FMI SAMPLE ID`, row.names(RNA_data))

# Preprocess DNA Data
DNA_data <- subset(DNA_data, select = c("FMI SAMPLE ID", "SAMPLE ID", "DISEASE"))
DNA_data <- DNA_data[DNA_data$`FMI SAMPLE ID` %in% common_fmis, ]
row.names(DNA_data) <- DNA_data$`FMI SAMPLE ID`
DNA_data <- DNA_data[common_fmis, ]
RNA_data <- RNA_data[common_fmis, ]
row.names(DNA_data) <- DNA_data$`FMI SAMPLE ID`

# Ensure matching row names
stopifnot(all(row.names(DNA_data) == row.names(RNA_data)))
row.names(DNA_data) <- DNA_data$`SAMPLE ID`
row.names(RNA_data) <- DNA_data$`SAMPLE ID`

DNA_data_UPS <- DNA_data[DNA_data$DISEASE == "Soft tissue sarcoma (NOS)",]

# Convert RNA data to numeric and transpose
RNA_data <- data.frame(lapply(RNA_data, as.numeric))
RNA_data <- data.frame(t(RNA_data))

# Filter rows with non-zero sums
RNA_data <- RNA_data[rowSums(RNA_data) > 0, ]
colnames(RNA_data) <- row.names(DNA_data)

ups_ids <- DNA_data_UPS$`SAMPLE ID`

RNA_data <- RNA_data[,ups_ids]

# Apply voom transformation
Voom <- voom(RNA_data, plot = FALSE, normalize.method = "quantile")

Voom <- data.frame(Voom$E)

clinical_data_clean <- read.csv("FILES/clinical_data_clean.csv", row.names=1)

# Remove outlier
Voom <- Voom[, colnames(Voom) != "X1.102.1"]
Voom <- data.frame(Voom)

RNA <- data.frame(voom(RNA_data, plot = FALSE,normalize.method = "quantile"))

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Sarcoma.Histopathological.Subtype"] <- "Sarcoma Histopathological Subtype"

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Neo.Adjuvant...Adjuvant.Treatment"] <- "Neo Adjuvant / Adjuvant Treatment"

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Local.Recurrence"] <- "Local Recurrence"

colnames(clinical_data_clean)[colnames(clinical_data_clean) == "Distant.Recurrence"] <- "Distant Recurrence"

colnames(RNA) <- sub("^X", "", colnames(RNA))

common_patients <- intersect(colnames(RNA),row.names(clinical_data_clean))

clinical_data_clean <- clinical_data_clean[common_patients,]
RNA <- RNA[,common_patients]

RNA <- data.frame(t(RNA))

# Prepare dataframe to store p-values
gene_p_values <- data.frame(Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)

# Perform survival analysis gene-by-gene
for (gene in colnames(RNA)) {
  
  print(gene)
  clinical_data_clean$gene_val <- as.numeric(RNA[, gene])
  
  # Categorize gene expression as "High" or "Low"
  mean_val <- mean(clinical_data_clean$gene_val, na.rm = TRUE)
  clinical_data_clean$gene_group <- ifelse(clinical_data_clean$gene_val > mean_val, "High", "Low")
  
  # Survival analysis (KM)
  surv_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY, Morte.S.N) ~ gene_group, data = clinical_data_clean)
  
  # Log-rank test
  p_val <- surv_pvalue(surv_fit)$pval
  
  # Extract final survival probabilities
  surv_summary <- summary(surv_fit)
  
  # Extract last survival probability of each group
  last_surv <- tapply(surv_summary$surv, surv_summary$strata, tail, 1)
  
  # Compare final probabilities directly
  high_worse <- last_surv["gene_group=High"] < last_surv["gene_group=Low"]
  
  # Append to dataframe
  gene_p_values <- rbind(gene_p_values, data.frame(Gene = gene,
                                                   P_Value = p_val,
                                                   High_Worse = high_worse))
  
  
  
}

rownames(gene_p_values) <- gene_p_values$Gene

write.csv(gene_p_values, "surv_p_values_per_gene_UPS.csv")


gene_p_values <- read_csv("RESULTS/surv_p_values_per_gene_UPS.csv")
gene_p_values_high <-gene_p_values[gene_p_values$High_Worse,]

gene_p_values_high <- gene_p_values_high[gene_p_values_high$P_Value<0.05,]

gene_p_values_high <- gene_p_values_high[!is.na(gene_p_values_high$P_Value), ]



library(org.Hs.eg.db)
go_enrich <- enrichGO(gene =   gene_p_values_high$Gene,
                      universe = colnames(RNA),
                      OrgDb = org.Hs.eg.db,
                      keyType="SYMBOL",
                      ont = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.05,
                      readable = TRUE, 
                      minGSSize = 3)

barplot(go_enrich, title = "Over Represented Pathways",showCategory = 20)
cnetplot(go_enrich,max.overlaps =200)


