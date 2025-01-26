library(readxl)

# Define output file path
output_file <- "FILES/clinical_data_clean.csv"

# Check if the output file exists
if (!file.exists(output_file)) {
  # Load clinical data
  clinical_data <- read.csv("FILES/Base de dados Clinica_update.csv", sep = ";")
  
  # Convert date columns
  clinical_data$Data.cirurgia <- as.Date(clinical_data$Data.cirurgia)
  clinical_data$Data.UF <- as.Date(clinical_data$Data.UF)
  
  # Keep only records where surgery data is not missing
  clinical_data <- clinical_data[!is.na(clinical_data$Data.cirurgia), ]
  
  # Calculate time in days from surgery to death/last follow-up
  clinical_data$TIME_DEATH_FROM_SURGERY <- as.numeric(clinical_data$Data.UF - clinical_data$Data.cirurgia)
  
  # Extract birth and diagnostic year
  clinical_data$Birth_Year <- substr(clinical_data$Data.nascimento, 1, 4)
  clinical_data$Diagnostic_Year <- substr(clinical_data$Data.diagnóstico., 1, 4)
  
  # Calculate age at the time of surgery
  clinical_data$Age <- as.numeric(clinical_data$Data.cirurgia) - as.numeric(clinical_data$Birth_Year)
  
  # Convert disease names from Portuguese to English
  clinical_data$Histologia <- ifelse(
    clinical_data$Histologia == "LIPOSSARCOMA DESDIFERENCIADO", "DDLPS",
    ifelse(clinical_data$Histologia == "LEIOMIOSSARCOMA", "LMS",
           ifelse(clinical_data$Histologia == "SARCOMA PLEIOMÓRFICO INDIFERENCIADO", "UPS",
                  "")
    )
  )
  
  # Rename columns to English
  colnames(clinical_data)[colnames(clinical_data) == "Histologia"] <- "Sarcoma Histopathological Subtype"
  colnames(clinical_data)[colnames(clinical_data) == "Tratamento.adjuvante"] <- "Neo Adjuvant / Adjuvant Treatment"
  colnames(clinical_data)[colnames(clinical_data) == "Sexo"] <- "Gender"
  colnames(clinical_data)[colnames(clinical_data) == "Recidiva.S.N"] <- "Local Recurrence"
  colnames(clinical_data)[colnames(clinical_data) == "Metastização.S.N"] <- "Distant Recurrence"
  
  
  # Set row names to study numbers
  row.names(clinical_data) <- clinical_data$Nº.Estudo
  
  # Remove the specific study
  clinical_data <- clinical_data[clinical_data$Nº.Estudo != "1.44.1", ]
  
  # Write the cleaned data to the output file
  write.csv(clinical_data, output_file)
  
  print("File has been processed and saved.")
} else {
  print("Output file already exists. No action taken.")
}
