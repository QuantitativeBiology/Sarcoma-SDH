library(EnhancedVolcano)
library(limma)
library(ComplexHeatmap)
library(gplots)
library(edgeR)
library(GSEABase)
library(fgsea)
#library(ssGSEA2)
library(corto)
library(survival)
library(survminer)
library(rempsyc)
library(msigdbr)
library(clusterProfiler)
library (enrichplot)
library(gplots)

library(rempsyc)

RNA_data <- read.csv("FILES/TCGA-SARC.csv",row.names=1)

clinical_data <- read.csv("FILES/TCGA-SARC-CLINICAL.csv",row.names=1)

table(clinical_data$PaperHistology)

unique(clinical_data$PaperHistology)

clinical_data <- clinical_data[clinical_data$PaperHistology %in% c("leiomyosarcoma - soft tissue","dedifferentiated liposarcoma", "undifferentiated pleomorphic sarcoma" ),]

clinical_data$Patient <- gsub("-", ".", clinical_data$Patient, fixed = TRUE)

RNA_data <- RNA_data[,colnames(RNA_data) %in% clinical_data$Patient]

clinical_data <- clinical_data[order(clinical_data$PaperHistology),]

RNA_data <- RNA_data[,clinical_data$Patient]


disease_ddlps <- as.numeric(grepl("leiomyosarcoma - soft tissue", clinical_data$PaperHistology, ignore.case = TRUE))
disease_lms <- as.numeric(grepl("dedifferentiated liposarcoma", clinical_data$PaperHistology, ignore.case = TRUE))
disease_ups <- as.numeric(grepl("undifferentiated pleomorphic sarcoma", clinical_data$PaperHistology, ignore.case = TRUE))
# 
# 
design <- cbind(disease_ddlps,disease_lms,disease_ups)
# 
keep <- filterByExpr(RNA_data, design = design)

RNA_data <- RNA_data[keep,]

#RNA_data <- DGEList(counts = RNA_data, genes = rownames(RNA_data))

# Normalize the counts using the TMM method
#RNA_data <- calcNormFactors(RNA_data, method = "TMM")

cont.matrix <- makeContrasts(disease_lms-disease_ddlps,
                             disease_ups-disease_lms,
                             disease_ups-disease_ddlps, levels=design)


Voom <- voom(RNA_data, design, plot = FALSE,normalize.method = "quantile")

voom_df <- data.frame(Voom)


vfit <- lmFit(Voom, design)
vfit  <- contrasts.fit(vfit,cont.matrix)
efit <- eBayes(vfit)

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
    write.csv(df_to_save, "RESULTS/deg_tcga_leiomyo_lipo.csv")
    deg_leiomyo_lipo <- df_to_save
    
    deg_leiomyo_lipo_degs <- df_to_save[df_to_save$adj.P.Val<0.05,]
    write.csv(deg_leiomyo_lipo_degs, "RESULTS/deg_leiomyo_lipo_sign_tcga.csv")
    
  } else if (coefi == 2) {
    title <- "UPS vs LMS"
    write.csv(df_to_save, "RESULTS/deg_tcga_ups_leiomyo.csv")
    deg_ups_leiomyo <- df_to_save
    deg_ups_leiomyo_sign <- df_to_save[df_to_save$adj.P.Val<0.05,]
    write.csv(deg_ups_leiomyo_sign, "RESULTS/deg_ups_leiomyo_sign_tcga.csv")
  } else if (coefi == 3) {
    title <- "UPS vs DDLPS"
    write.csv(df_to_save, "RESULTS/deg_tcga_ups_lipo.csv")
    deg_ups_lipo <- df_to_save
    deg_ups_lipo_degs <- df_to_save[df_to_save$adj.P.Val<0.05,]
    write.csv(deg_ups_lipo_degs, "RESULTS/deg_ups_ddlps_sign_tcga.csv")
  }
  print(title)
  for (i in 1:length(rownames(deg))) {
    if (y[i] <= 0.05 & x[i] %in% metabolic_genes$gene_symbol) {
      #x[i] <- ""
      df <- rbind(df, c(x[i],y[i], log[i]))
    } else{
      x[i] <- ""
    }
  }
  colnames(df) <- names
  row.names(df) <- df$RowNames
  
  vulcano <- EnhancedVolcano(
    deg,
    lab = x,
    x = 'logFC',
    y = 'P.Value',
    xlim = c(-8,8),
    ylim = c(0,30),
    labSize = 3.0,
    FCcutoff = 0.85,
    pCutoff = 1e-02,
    title = title,
    legendPosition = 'none',
    drawConnectors = T,
    max.overlaps = 100
  )
  vulcano_plots[[coefi]] <- vulcano

  sign <- df_to_save[df_to_save$adj.P.Val <= 0.05,]

  print("SDHB" %in% row.names(sign))

  deg <- deg[deg$adj.P.Val<0.05,]

  print("Mean of LOGFC")
  print(mean(deg$logFC))

  deg_positives <- deg[deg$logFC>0.0,]
  deg_negatives <- deg[deg$logFC<0.0,]
  #print(deg_positives)

  print("Number of DEGS")
  print(length(row.names(deg_positives)) + length(row.names(deg_negatives)))

  print("Mean of over")
  print(mean(deg_positives$logFC))

  print("Mean of under")
  print(mean(deg_negatives$logFC))
}
x <- vulcano_plots[[1]]
y <- vulcano_plots[[2]]
z <- vulcano_plots[[3]]
gridExtra::grid.arrange(x,y,z, ncol = 3)

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
  
  ups_lms_p_value <- as.character(round(ups_leiomyo$adj.P.Val, 3))
  if ( ups_lms_p_value == 0){
    ups_lms_p_value <- "<0.001"
  }
  
  ups_ddlps_p_value <- as.character(round(ups_lipo$adj.P.Val, 3))
  if ( ups_ddlps_p_value == 0){
    ups_ddlps_p_value <- "<0.001"
  }

  # Extract and round adj_p_values and logfc_values
  adj_p_values <- c(
    as.character(round(leiomyo_lipo$adj.P.Val, 3)), 
    ups_lms_p_value, 
    ups_ddlps_p_value
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

all(colnames(voom_df) == clinical_data$Patient)

voom_df <- data.frame(t(voom_df))

voom_df$subtype <- clinical_data$PaperHistology


voom_df$subtype <- ifelse(voom_df$subtype  == "dedifferentiated liposarcoma", "DDLPS", voom_df$subtype)
voom_df$subtype <- ifelse(voom_df$subtype  == "leiomyosarcoma - soft tissue" , "LMS", voom_df$subtype)
voom_df$subtype <- ifelse(voom_df$subtype  == "undifferentiated pleomorphic sarcoma" , "UPS", voom_df$subtype)

voom_long <- voom_df %>%
  dplyr::select(subtype, SDHA,SDHB, SDHC, SDHD) %>%
  tidyr::pivot_longer(
    cols = c(SDHA,SDHB, SDHC, SDHD),   # The genes you want to plot
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

voom_df$TIME_DEATH_FROM_SURGERY <- clinical_data$Last_FU

voom_df$DEATH <- clinical_data$Status

voom_df$original_class <- clinical_data$Histology

voom_df$Paper_Histology <- clinical_data$PaperHistology

voom_df$FNCLCC_GRADE <- as.character(clinical_data$FNCLCC_GRADE)

voom_df$SDHB_exp <- voom_df$SDHB

voom_df$Gender <- clinical_data$Gender

# All Subtypes
voom_df$SDHB <- ifelse(voom_df$SDHB_exp > mean(voom_df$SDHB_exp), "High", "Low")

km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY,DEATH ) ~ SDHB, data=voom_df)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival SDHB")

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, DEATH) ~ FNCLCC_GRADE + Paper_Histology + SDHB + Gender , data = voom_df)
ggforest(cox, fontsize = 1)

# Just UPS
voom_df <- voom_df[voom_df$subtype == "UPS",]

voom_df$SDHB <- ifelse(voom_df$SDHB_exp > mean(voom_df$SDHB_exp), "High", "Low")

km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY,DEATH) ~ SDHB, data=voom_df)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival SDHB: Only UPS")
x

cox <- coxph(Surv(TIME_DEATH_FROM_SURGERY, DEATH) ~ FNCLCC_GRADE + SDHB + Gender , data = voom_df)
ggforest(cox, fontsize = 1)


