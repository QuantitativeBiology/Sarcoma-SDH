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
  } else if (coefi == 2) {
    title <- "UPS vs LMS"
    write.csv(df_to_save, "RESULTS/deg_tcga_ups_leiomyo.csv")
  } else if (coefi == 3) {
    title <- "UPS vs DDLPS"
    write.csv(df_to_save, "RESULTS/deg_tcga_ups_lipo.csv")
  }
  print(title)
  for (i in 1:length(rownames(deg))) {
    if (y[i] <= 0.05) {
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
    drawConnectors = F,
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


all(colnames(voom_df) == clinical_data$Patient)

voom_df <- data.frame(t(voom_df))

voom_df$subtype <- clinical_data$PaperHistology


voom_df$subtype <- ifelse(voom_df$subtype  == "dedifferentiated liposarcoma", "DDLPS", voom_df$subtype)
voom_df$subtype <- ifelse(voom_df$subtype  == "leiomyosarcoma - soft tissue" , "LMS", voom_df$subtype)
voom_df$subtype <- ifelse(voom_df$subtype  == "undifferentiated pleomorphic sarcoma" , "UPS", voom_df$subtype)

voom_long <- voom_df %>%
  dplyr::select(subtype, SDHB, SDHC, SDHD) %>%
  tidyr::pivot_longer(
    cols = c(SDHB, SDHC, SDHD),   # The genes you want to plot
    names_to = "Gene",            # New column that will store gene names
    values_to = "Expression"      # New column that will store expression values
  )


p <- ggplot(voom_long, aes(x = subtype, y = Expression, fill = subtype)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  theme_pubr() +
  labs(x = "", y = "Expression (Voom Normalized)") +
  facet_wrap(~Gene, ncol = 3) +            # Facet by Gene; 3 columns side by side
  theme(legend.position = "none")          # Remove legend if not needed

print(p)

voom_df$TIME_DEATH_FROM_SURGERY <- clinical_data$Last_FU

voom_df$DEATH <- clinical_data$Status

voom_df$original_class <- clinical_data$Histology

voom_df$SDHB_exp <- voom_df$SDHB

# All Subtypes
voom_df$SDHB <- ifelse(voom_df$SDHB_exp > mean(voom_df$SDHB_exp), "High", "Low")

km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY,DEATH ) ~ SDHB, data=voom_df)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival SDHB")

x

# Just UPS
voom_df <- voom_df[voom_df$subtype == "UPS",]

voom_df$SDHB <- ifelse(voom_df$SDHB_exp > mean(voom_df$SDHB_exp), "High", "Low")

km_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY,DEATH) ~ SDHB, data=voom_df)

x <- ggsurvplot(km_fit,pval=TRUE,risk.table=TRUE, conf.int = TRUE, 
                title = "Overall Survival SDHB")
x

# Count the number of "High" occurrences for each original class
sdhb_counts <- table(voom_df$original_class[voom_df$SDHB == "High"])

# Convert to dataframe for ggplot
sdhb_df <- as.data.frame(sdhb_counts)
colnames(sdhb_df) <- c("Original_Class", "High_Count")

# Plot using ggplot2
ggplot(sdhb_df, aes(x = reorder(Original_Class, -High_Count), y = High_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Number of SDHB High Cases per Original Class",
       x = "Original Class",
       y = "Count of SDHB High") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Count the number of "High" occurrences for each original class
sdhb_counts <- table(voom_df$original_class[voom_df$SDHB == "Low"])

# Convert to dataframe for ggplot
sdhb_df <- as.data.frame(sdhb_counts)
colnames(sdhb_df) <- c("Original_Class", "High_Count")

# Plot using ggplot2
ggplot(sdhb_df, aes(x = reorder(Original_Class, -High_Count), y = High_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Number of SDHB Low Cases per Original Class",
       x = "Original Class",
       y = "Count of SDHB Low") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
