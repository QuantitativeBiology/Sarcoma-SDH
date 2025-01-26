# Main script to execute all steps sequentially

# Define the list of scripts in execution order
scripts <- c(
  "0.clinical_data_preprocessing.R",
  "0.rna_normalization.R"
)

# Sequentially source each script
for (script in scripts) {
  cat("Running script:", script, "\n")
  source(script)
  cat("Finished running script:", script, "\n\n")
}

cat("All scripts executed successfully.\n")