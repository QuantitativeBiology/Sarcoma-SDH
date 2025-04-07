# Load necessary libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# Query gene expression data
query <- GDCquery(project = "TCGA-SARC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  sample.type = "Primary Tumor")

# Download data
GDCdownload(query)

# Prepare as SummarizedExperiment object
sarcoma_se <- GDCprepare(query, summarizedExperiment = TRUE)

# Check metadata
colnames(colData(sarcoma_se))

data <- data.frame(sarcoma_se)

ups <- sarcoma_se[sarcoma_se$paper_histology == "undifferentiated pleomorphic sarcoma",]

# Extract clinical metadata from the object
clinical <- data.frame(
  Patient = colData(sarcoma_se)$barcode,
  Histology = colData(sarcoma_se)$primary_diagnosis,
  Last_FU = colData(sarcoma_se)$days_to_last_follow_up,
  Status = colData(sarcoma_se)$vital_status,
  Tissue_Type = colData(sarcoma_se)$tissue_or_organ_of_origin,
  Grade = colData(sarcoma_se)$tumor_grade,
  FNCLCC_Grade = colData(sarcoma_se)$paper_FNCLCC_grade,
  Gender = colData(sarcoma_se)$gender,
  Progression_Or_Recurrence = colData(sarcoma_se)$progression_or_recurrence,
  PaperHistology = colData(sarcoma_se)$paper_histologic_diagnosis,
  SiteOfResection = colData(sarcoma_se)$site_of_resection_or_biopsy
)

# Check clinical info
head(clinical)
