library(clusterProfiler)
library(ggplot2)
library(plotly)
library(org.Hs.eg.db)

dge_ups_lipo <- read.csv("RESULTS/deg_ups_lipo.csv", row.names=1)
dge_ups_lipo$genes <- NULL
names <- row.names(dge_ups_lipo)
dge_ups_lipo <- data.frame(sapply(dge_ups_lipo, as.numeric))
row.names(dge_ups_lipo) <- names


dge_ups_lipo_sign <- dge_ups_lipo[dge_ups_lipo$adj.P.Val <= 0.05,]
dge_ups_lipo_over <- dge_ups_lipo_sign[dge_ups_lipo_sign$logFC > 0,]
dge_ups_lipo_under <- dge_ups_lipo_sign[dge_ups_lipo_sign$logFC< 0,]

dge_ups_leiomyo <- read.csv("RESULTS/deg_ups_leiomyo.csv", row.names = 1)
dge_ups_leiomyo$genes <- NULL
names <- row.names(dge_ups_leiomyo)
dge_ups_leiomyo <- data.frame(sapply(dge_ups_leiomyo, as.numeric))
row.names(dge_ups_leiomyo) <- names

dge_ups_leiomyo_sign <- dge_ups_leiomyo[dge_ups_leiomyo$adj.P.Val <= 0.05,]
dge_ups_leiomyo_over <- dge_ups_leiomyo_sign[dge_ups_leiomyo_sign$logFC > 0,]
dge_ups_leiomyo_under <- dge_ups_leiomyo_sign[dge_ups_leiomyo_sign$logFC < 0,]

dge_leiomyo_lipo <- read.csv("RESULTS/deg_leiomyo_lipo.csv", row.names = 1)
dge_leiomyo_lipo$genes <- NULL
names <- row.names(dge_leiomyo_lipo)
dge_leiomyo_lipo <- data.frame(sapply(dge_leiomyo_lipo, as.numeric))
row.names(dge_leiomyo_lipo) <- names
dge_leiomyo_lipo_sign <- dge_leiomyo_lipo[dge_leiomyo_lipo$adj.P.Val <= 0.05,]
dge_leiomyo_lipo_over <- dge_leiomyo_lipo_sign[dge_leiomyo_lipo_sign$logFC > 0,]
dge_leiomyo_lipo_under <- dge_leiomyo_lipo_sign[dge_leiomyo_lipo_sign$logFC < 0,]

go_enrich_ups_lipo_over <- enrichGO(gene = row.names(dge_ups_lipo_over),
                                    universe = row.names(dge_ups_lipo),
                                    OrgDb = org.Hs.eg.db,
                                    keyType="SYMBOL",
                                    ont = "BP",
                                    pAdjustMethod = "fdr",
                                    pvalueCutoff = 0.05,
                                    minGSSize = 3,
                                    readable = TRUE)


go_enrich_ups_leiomyo_over <- enrichGO(gene = row.names(dge_ups_leiomyo_over),
                                       universe = row.names(dge_ups_leiomyo),
                                       OrgDb = org.Hs.eg.db,
                                       keyType="SYMBOL",
                                       ont = "BP",
                                       pAdjustMethod = "fdr",
                                       pvalueCutoff = 0.05,
                                       minGSSize = 3,
                                       readable = TRUE)

go_enrich_leiomyo_lipo_over <-  enrichGO(gene = row.names(dge_leiomyo_lipo_over),
                                           universe = row.names(dge_leiomyo_lipo),
                                           OrgDb = org.Hs.eg.db,
                                           keyType="SYMBOL",
                                           ont = "BP",
                                           pAdjustMethod = "fdr",
                                           pvalueCutoff = 0.05,
                                           minGSSize = 3,
                                           readable = TRUE)

barplot(go_enrich_ups_lipo_over, title= "UPS vs DDLPS", showCategory = 20)
barplot(go_enrich_ups_leiomyo_over, title= "UPS vs LMS", showCategory = 20)
barplot(go_enrich_leiomyo_lipo_over, title= "LMS vs DDLPS", showCategory = 20)

# Get metabolic genes and pathways 
metabolic_genes <- read.csv("RESULTS/METABOLIC_GENES_IN_KIT.csv")

#metabolic_genes <- metabolic_genes[metabolic_genes$gene_symbol %in% c("SDHA","SDHB", "SDHC", "SDHD"),]

x <- rownames(dge_leiomyo_lipo)
y <- dge_leiomyo_lipo$adj.P.Val
df <- data.frame(RowNames = character(), Value = numeric())
names <- colnames(df)

for (i in 1:length(rownames(dge_leiomyo_lipo))) {
  if (!is.na(y[i]) && y[i] <= 0.05 && x[i] %in% metabolic_genes$gene_symbol) {
    #x[i] <- ""
    df <- rbind(df, c(x[i],y[i]))
  } else{
    x[i] <- ""
  }
}

leiomyo_lipo <- EnhancedVolcano(
  dge_leiomyo_lipo,
  lab = x,
  x = 'logFC',
  y = 'P.Value',
  xlim = c(-8,8),
  ylim = c(0,30),
  labSize = 3.0,
  FCcutoff = 0.75,
  pCutoff = 1e-02,
  title = 'LMS vs LMS',
  legendPosition = 'none',
  drawConnectors = TRUE,
  max.overlaps = 100
)

x <- rownames(dge_ups_leiomyo)
y <- dge_ups_leiomyo$adj.P.Val
df <- data.frame(RowNames = character(), Value = numeric())
names <- colnames(df)

for (i in 1:length(rownames(dge_ups_leiomyo))) {
  if (!is.na(y[i]) && y[i] <= 0.05 && x[i] %in% metabolic_genes$gene_symbol) {
    #x[i] <- ""
    df <- rbind(df, c(x[i],y[i]))
  } else{
    x[i] <- ""
  }
}

ups_leiomyo <- EnhancedVolcano(
  dge_ups_leiomyo,
  lab = x,
  x = 'logFC',
  y = 'P.Value',
  xlim = c(-8,8),
  ylim = c(0,30),
  labSize = 3.0,
  FCcutoff = 0.75,
  pCutoff = 1e-02,
  title = 'UPS vs LMS',
  legendPosition = 'none',
  drawConnectors = TRUE,
  max.overlaps = 100
)


x <- rownames(dge_ups_lipo)
y <- dge_ups_lipo$adj.P.Val
df <- data.frame(RowNames = character(), Value = numeric())
names <- colnames(df)

for (i in 1:length(rownames(dge_ups_lipo))) {
  if (!is.na(y[i]) && y[i] <= 0.05 && x[i] %in% metabolic_genes$gene_symbol) {
    #x[i] <- ""
    df <- rbind(df, c(x[i],y[i]))
  } else{
    x[i] <- ""
  }
}

ups_lipo <- EnhancedVolcano(
  dge_ups_lipo,
  lab = x,
  x = 'logFC',
  y = 'P.Value',
  xlim = c(-8,8),
  ylim = c(0,30),
  labSize = 3.0,
  FCcutoff = 0.75,
  pCutoff = 1e-02,
  title = 'UPS vs DDLPS',
  legendPosition = 'none',
  drawConnectors = TRUE,
  max.overlaps = 100
)

leiomyo_lipo

ups_leiomyo

ups_lipo


gridExtra::grid.arrange(leiomyo_lipo,ups_leiomyo,ups_lipo, ncol = 3)






