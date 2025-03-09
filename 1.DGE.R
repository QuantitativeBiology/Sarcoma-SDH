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

# Define output file path
output_file <- "FILES/normalized_counts.csv"


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

# )
cont.matrix <- makeContrasts(disease_leimyo-disease_lipo,
                             disease_nos-disease_leimyo,
                             disease_nos-disease_lipo, levels=design)

Voom <- voom(RNA_data, plot = FALSE,normalize.method = "quantile")

vfit <- lmFit(Voom, design)
vfit  <- contrasts.fit(vfit,cont.matrix)
efit <- eBayes(vfit)
plotSA(efit, main = 'final model: Mean-Variance trend')

metabolic_genes <- read.csv("RESULTS/METABOLIC_GENES_IN_KIT.csv")

deg_leiomyo_lipo <- NULL
deg_ups_leiomyo <- NULL
deg_ups_lipo <- NULL

vulcano_plots <- list()
for(coefi in 1:3){
  deg <- topTable(efit, coef = coefi,adjust.method = 'fdr', number=Inf)
  x <- rownames(deg)
  y <- deg$adj.P.Val
  log <- deg$logFC
  df <- data.frame(RowNames = character(), Adj_P = numeric(), LogFc = numeric())
  df_to_save <- data.frame(deg)
  names <- colnames(df)
  
  title = ''
  if (coefi == 1) {
    title <- "LMS vs DDLPS"
    deg_leiomyo_lipo <- df_to_save
    write.csv(df_to_save, "RESULTS/deg_leiomyo_lipo.csv")
    deg_leiomyo_lipo_degs <- df_to_save[df_to_save$adj.P.Val<0.05,]
    write.csv(deg_leiomyo_lipo_degs, "RESULTS/deg_leiomyo_lipo_sign.csv")
  } else if (coefi == 2) {
    title <- "UPS vs LMS"
    deg_ups_leiomyo <- df_to_save
    write.csv(df_to_save, "RESULTS/deg_ups_leiomyo.csv")
    deg_ups_leiomyo_sign <- df_to_save[df_to_save$adj.P.Val<0.05,]
    write.csv(deg_ups_leiomyo_sign, "RESULTS/deg_ups_leiomyo_sign.csv")
  } else if (coefi == 3) {
    title <- "UPS vs DDLPS"
    deg_ups_lipo <- df_to_save
    deg_ups_lipo_degs <- df_to_save[df_to_save$adj.P.Val<0.05,]
    write.csv(deg_ups_lipo_degs, "RESULTS/deg_ups_ddlps_sign.csv")
  }
  for (i in 1:length(rownames(deg))) {
    if (y[i] <= 0.05){ #&& x[i] %in% metabolic_genes$gene_symbol
      #x[i] <- ""
      df <- rbind(df, c(x[i],y[i], log[i]))
    } else{
      x[i] <- ""
    }
  }
  colnames(df) <- names
  row.names(df) <- df$RowNames
  print(title)
  print(length(row.names(df)))
  
  print(mean(deg$logFC))
  
  vulcano <- EnhancedVolcano(
    deg,
    lab = x,
    x = 'logFC',
    y = 'P.Value',
    xlim = c(-8,8),
    ylim = c(0,30),
    labSize = 3.0,
    FCcutoff = 0.5,
    pCutoff = 1e-02,
    title = title,
    legendPosition = 'none',
    drawConnectors = T,
    max.overlaps = 100
  )
  vulcano_plots[[coefi]] <- vulcano
  
  deg <- deg[deg$adj.P.Val<0.05,]
  
  deg_positives <- deg[deg$logFC>0.0,]
  deg_negatives <- deg[deg$logFC<0.0,]
  
  print(mean(deg_positives$logFC))
  print(mean(deg_negatives$logFC))
  
}
x <- vulcano_plots[[1]]
y <- vulcano_plots[[2]]
z <- vulcano_plots[[3]]
gridExtra::grid.arrange(x,y,z, ncol = 3)

voom_dataframe <- data.frame(Voom)
colnames(voom_dataframe) <- sub("^X", "", colnames(voom_dataframe))

voom_dataframe <- data.frame(t(voom_dataframe))

DNA_data$DISEASE <- ifelse(DNA_data$DISEASE == "Soft tissue sarcoma (NOS)", "UPS", DNA_data$DISEASE)
DNA_data$DISEASE <- ifelse(DNA_data$DISEASE == "Soft tissue leiomyosarcoma", "LMS", DNA_data$DISEASE)
DNA_data$DISEASE <- ifelse(DNA_data$DISEASE == "Soft tissue liposarcoma", "DDLPS", DNA_data$DISEASE)

voom_dataframe$subtype <- DNA_data$DISEASE

# Ensure subtype is treated as a factor (if it's not already)
voom_dataframe$subtype <- factor(voom_dataframe$subtype)

# Create the violin plot
p <- ggplot(voom_dataframe, aes(x = subtype, y = SDHB, fill = subtype)) +
  geom_violin(trim = FALSE, alpha = 0.7) +       # Violin layer
  geom_boxplot(width = 0.2, fill = "white",      # Boxplot on top (optional)
               outlier.shape = NA) +
  # Optionally, add individual points or jitter:
  # geom_jitter(width = 0.15, size = 1, alpha = 0.5)
  theme_pubr() +                                 # Publication-like theme from ggpubr
  labs(x = "", y = "SDHB Expression (Voom Normalized)") +    # Axis labels
  theme(legend.position = "none")                # Remove legend if not needed

# Display the plot
print(p)


# Create the violin plot
p <- ggplot(voom_dataframe, aes(x = subtype, y = SDHC, fill = subtype)) +
  geom_violin(trim = FALSE, alpha = 0.7) +       # Violin layer
  geom_boxplot(width = 0.2, fill = "white",      # Boxplot on top (optional)
               outlier.shape = NA) +
  # Optionally, add individual points or jitter:
  # geom_jitter(width = 0.15, size = 1, alpha = 0.5)
  theme_pubr() +                                 # Publication-like theme from ggpubr
  labs(x = "", y = "SDHC Expression (Voom Normalized)") +    # Axis labels
  theme(legend.position = "none")                # Remove legend if not needed

# Display the plot
print(p)

p <- ggplot(voom_dataframe, aes(x = subtype, y = SDHD, fill = subtype)) +
  geom_violin(trim = FALSE, alpha = 0.7) +       # Violin layer
  geom_boxplot(width = 0.2, fill = "white",      # Boxplot on top (optional)
               outlier.shape = NA) +
  # Optionally, add individual points or jitter:
  # geom_jitter(width = 0.15, size = 1, alpha = 0.5)
  theme_pubr() +                                 # Publication-like theme from ggpubr
  labs(x = "", y = "SDHD Expression (Voom Normalized)") +    # Axis labels
  theme(legend.position = "none")                # Remove legend if not needed

# Display the plot
print(p)

# Convert voom_dataframe from wide to long format
voom_long <- voom_dataframe %>%
  dplyr::select(subtype, SDHA, SDHB, SDHC, SDHD) %>%
  tidyr::pivot_longer(
    cols = c(SDHA, SDHB, SDHC, SDHD),   # The genes you want to plot
    names_to = "Gene",            # New column that will store gene names
    values_to = "Expression"      # New column that will store expression values
  )


# Define a colorblind-friendly palette
palette <- brewer.pal(8, "Set2")

# Enhanced plot for Nature publication with fixed y-axis
p <- ggplot(voom_long, aes(x = subtype, y = Expression, fill = subtype)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black", size = 0.4) +  # Adjust alpha and outline
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, color = "black", size = 0.3) +
  scale_fill_manual(values = palette) +  # Apply custom color palette
  theme_classic() +  # Use a clean theme
  labs(x = "", y = "Expression (Voom Normalized)") +
  facet_wrap(~Gene, ncol = 4, scales = "fixed") +  # Fixed y-axis across plots
  theme(
    legend.position = "none",  # Remove legend if not needed
    strip.text = element_text(size = 10, face = "bold"),  # Facet labels
    axis.text = element_text(size = 9),  # Axis labels
    axis.title = element_text(size = 10, face = "bold"),  # Axis titles
    plot.title = element_text(size = 12, face = "bold"),  # Plot title
    axis.line = element_line(size = 0.5, color = "black"),  # Axis lines
    strip.background = element_blank(),  # Clean facet background
    panel.spacing = unit(0.5, "lines")  # Adjust spacing between facets
  )

print(p)

# List of genes to extract
genes <- c("SDHA", "SDHB", "SDHC", "SDHD")

# Prepare a list to store dataframes for each gene
deg_tables <- list()

# Iterate over each gene and extract data
for (gene in genes) {
  # Extract rows for the current gene
  deg_leiomyo_lipo$genes <- row.names(deg_leiomyo_lipo)
  deg_ups_leiomyo$genes <- row.names(deg_ups_leiomyo)
  deg_ups_lipo$genes <- row.names(deg_ups_lipo)
  
  leiomyo_lipo <- deg_leiomyo_lipo[deg_leiomyo_lipo$genes == gene,]
  ups_leiomyo <- deg_ups_leiomyo[deg_ups_leiomyo$genes == gene,]
  ups_lipo <- deg_ups_lipo[deg_ups_lipo$genes == gene,]
  
  # Extract and round adj_p_values and logfc_values
  adj_p_values <- c(
    as.character(round(leiomyo_lipo$adj.P.Val, 3)), 
    as.character(round(ups_leiomyo$adj.P.Val, 3)), 
    as.character(round(ups_lipo$adj.P.Val, 3))
  )
  
  logfc_values <- c(
    as.character(round(leiomyo_lipo$logFC, 3)), 
    as.character(round(ups_leiomyo$logFC, 3)), 
    as.character(round(ups_lipo$logFC, 3))
  )
  
  # Create a dataframe for the current gene
  df <- data.frame(
    Comparison = c("LMS VS DDLPS", "UPS vs LMS", "UPS vs DDLPS"),
    adj_p_value = adj_p_values,
    logfc = logfc_values
  )
  
  # Store the dataframe in the list
  deg_tables[[gene]] <- df
  
  # Print the table using nice_table
  print(paste(gene, "Differential Gene Expression"))
  x <- nice_table(
    df,
    title = paste(gene, "Differential Gene Expression"),
    note = c("Limma-Voom Differential expression results")
  )
  
  print(x)
}
