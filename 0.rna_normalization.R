# Load Libraries
library(dplyr)
library(readxl)
library(limma)
library(edgeR)

# Define output file path
output_file <- "FILES/normalized_counts.csv"

# Check if the output file exists
if (!file.exists(output_file)) {
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
  
  # Convert RNA data to numeric and transpose
  RNA_data <- data.frame(lapply(RNA_data, as.numeric))
  RNA_data <- data.frame(t(RNA_data))
  
  # Filter rows with non-zero sums
  RNA_data <- RNA_data[rowSums(RNA_data) > 0, ]
  colnames(RNA_data) <- row.names(DNA_data)
  
  # Create binary columns for disease subtypes
  disease_nos <- as.numeric(grepl("NOS", DNA_data$DISEASE, ignore.case = TRUE))
  disease_lipo <- as.numeric(grepl("Soft tissue liposarcoma", DNA_data$DISEASE, ignore.case = TRUE))
  disease_leimyo <- as.numeric(grepl("Soft tissue leiomyosarcoma", DNA_data$DISEASE, ignore.case = TRUE))
  
  # Combine into design matrix
  design <- cbind(disease_lipo, disease_leimyo, disease_nos)
  
  # Filter RNA data by expression
  keep <- filterByExpr(RNA_data, design = design)
  RNA_data <- RNA_data[keep, ]
  
  # Apply voom transformation
  Voom <- voom(RNA_data, plot = FALSE, normalize.method = "quantile")
  
  Voom <- data.frame(Voom$E)
  
  # Remove outlier
  Voom <- Voom[, colnames(Voom) != "X1.102.1"]
  Voom <- data.frame(Voom)
  
  # Write the results to a CSV file
  write.csv(Voom, output_file)
  
  print("File has been processed and saved.")
} else {
  print("Output file already exists. No action taken.")
}
