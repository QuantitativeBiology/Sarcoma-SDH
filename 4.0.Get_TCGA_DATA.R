# Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(org.Hs.eg.db)

# Query TCGA-SARC RNA-Seq data
query_TCGA <- GDCquery(
  project = "TCGA-SARC",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  access = "open"
)

# Get results and download data
res <- getResults(query_TCGA)
GDCdownload(query_TCGA)

# Prepare SummarizedExperiment object
tcga_sarc_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)

# Obtain count matrix
sarc_matrix <- assay(tcga_sarc_data, "unstranded")
sarc_matrix <- data.frame(sarc_matrix)

# Extract clinical data from colData
clinical_data <- colData(tcga_sarc_data)


clinical_data



clinical_data$META=ifelse(clinical_data$`paper_distant recurrence`=='Distant Metastasis',1,0)
clinical_data$META[is.na(clinical_data$META)]=0
clinical_data$META_DATE=clinical_data$`paper_OS days`
clinical_data$META_DATE[which(clinical_data$META==1)]=clinical_data$`paper_distant recurrence days`[which(clinical_data$META==1)]





clinical <- data.frame(
  barcode                  = clinical_data$barcode,
  patient                  = clinical_data$patient,
  days_to_last_follow_up  = clinical_data$`paper_OS days`,
  days_to_death            = clinical_data$days_to_death,
  vital_status             = clinical_data$vital_status,
  sample_type              = clinical_data$sample_type,
  classification_of_tumor = clinical_data$classification_of_tumor,
  tissue_type              = clinical_data$tissue_type,
  size_tumor               = clinical_data$`paper_pathologic tumor size`,
  focality                 = clinical_data$`paper_disease multifocal indicator`,
  age_diagnosis            = clinical_data$`paper_age at diagnosis`,
  site_of_resection        = clinical_data$site_of_resection_or_biopsy,
  distant_recurrence       = clinical_data$META,
  distant_recurrence_time  = clinical_data$META_DATE
)

# Export clinical data to CSV
write.csv(clinical, "FILES/TCGA_clinical_data_updated.csv", row.names = FALSE)
