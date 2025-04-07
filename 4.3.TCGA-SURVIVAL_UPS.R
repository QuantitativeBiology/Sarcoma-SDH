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

clinical_data <- clinical_data[clinical_data$PaperHistology %in% c("undifferentiated pleomorphic sarcoma" ),]

clinical_data$Patient <- gsub("-", ".", clinical_data$Patient, fixed = TRUE)

RNA_data <- RNA_data[,colnames(RNA_data) %in% clinical_data$Patient]

clinical_data <- clinical_data[order(clinical_data$PaperHistology),]

RNA_data <- RNA_data[,clinical_data$Patient]

disease_ups <- as.numeric(grepl("undifferentiated pleomorphic sarcoma", clinical_data$PaperHistology, ignore.case = TRUE))


Voom <- voom(RNA_data, plot = FALSE,normalize.method = "quantile")

voom_df <- data.frame(Voom)

all(colnames(voom_df) == clinical_data$Patient)

voom_df <- data.frame(t(voom_df))

voom_df$subtype <- clinical_data$PaperHistology

voom_df$subtype <- ifelse(voom_df$subtype  == "undifferentiated pleomorphic sarcoma" , "UPS", voom_df$subtype)

voom_df$TIME_DEATH_FROM_SURGERY <- clinical_data$Last_FU

voom_df$DEATH <- clinical_data$Status

gene_p_values <- data.frame(Gene = character(), P_Value = numeric(), High_Worse = logical(), stringsAsFactors = FALSE)

#voom_df <- voom_df[,c("SDHB", "SDHC", "SDHD","TIME_DEATH_FROM_SURGERY", "DEATH", "subtype")]

# Loop gene-by-gene
for (gene in colnames(voom_df)) {
  if (!(gene %in% c("TIME_DEATH_FROM_SURGERY", "DEATH", "subtype"))){
    print(gene)
    voom_df$gene_val <- as.numeric(voom_df[, gene])
    
    # Categorize expression
    mean_val <- mean(voom_df$gene_val, na.rm = TRUE)
    voom_df$gene_group <- ifelse(voom_df$gene_val > mean_val, "High", "Low")
    
    # Check both groups exist
    if (length(unique(voom_df$gene_group)) < 2) next
    
    # Survival fit
    surv_fit <- survfit(Surv(TIME_DEATH_FROM_SURGERY, DEATH) ~ gene_group, data = voom_df)
    
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
}

gene_p_values <- read_csv("RESULTS/surv_p_values_per_gene_UPS.csv")

gene_p_values_high <-gene_p_values[!gene_p_values$High_Worse,]

gene_p_values_high <- gene_p_values_high[gene_p_values_high$P_Value<0.05,]

gene_p_values_high <- gene_p_values_high[!is.na(gene_p_values_high$P_Value), ]


gene_p_values_high$Gene

colnames(voom_df)


library(org.Hs.eg.db)
go_enrich <- enrichGO(gene =   gene_p_values_high$Gene,
                      universe = colnames(voom_df),
                      OrgDb = org.Hs.eg.db,
                      keyType="SYMBOL",
                      ont = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.05,
                      readable = TRUE, 
                      minGSSize = 10)

barplot(go_enrich, title = "Over Represented Pathways",showCategory = 20)
cnetplot(go_enrich,max.overlaps =200)






