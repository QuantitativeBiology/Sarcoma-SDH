library(clusterProfiler)
library(ggplot2)
library(plotly)

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

go_enrich_nos_lipo_over <- enrichGO(gene = row.names(dge_ups_lipo_over),
                                    universe = row.names(dge_ups_lipo),
                                    OrgDb = org.Hs.eg.db,
                                    keyType="SYMBOL",
                                    ont = "BP",
                                    pAdjustMethod = "fdr",
                                    pvalueCutoff = 0.05,
                                    minGSSize = 3,
                                    readable = TRUE)


go_enrich_nos_leiomyo_over <- enrichGO(gene = row.names(dge_ups_leiomyo_over),
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

barplot(go_enrich_nos_lipo_over, title= "UPS vs DDLPS", showCategory = 20)
barplot(go_enrich_nos_leiomyo_over, title= "UPS vs LMS", showCategory = 20)
barplot(go_enrich_leiomyo_lipo_over, title= "LMS vs DDLPS", showCategory = 20)

