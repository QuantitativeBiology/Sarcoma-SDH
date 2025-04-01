# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
library(org.Hs.eg.db)

# Load datasets
deg_ups_leiomyo <- read_csv("RESULTS/deg_ups_leiomyo.csv")
deg_ups_lipo <- read_csv("RESULTS/deg_ups_lipo.csv")
deg_leiomyo_lipo <- read_csv("RESULTS/deg_leiomyo_lipo.csv")

deg_ups_lipo_sign <- deg_ups_lipo[deg_ups_lipo$adj.P.Val<0.05,]
deg_ups_lipo_sign <- deg_ups_lipo_sign[deg_ups_lipo_sign$logFC>0.0,]
go_enrich1 <- enrichGO(gene =   deg_ups_lipo_sign$...1,
                      universe = deg_ups_lipo$...1,
                      OrgDb = org.Hs.eg.db,
                      keyType="SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.05,
                      readable = TRUE, 
                      minGSSize = 3)

barplot(go_enrich, title = "Over Represented Pathways UPS vs DDLPS",showCategory = 20)
cnetplot(go_enrich,max.overlaps =200)


deg_ups_leiomyo_sign <- deg_ups_leiomyo[deg_ups_leiomyo$adj.P.Val<0.05,]
deg_ups_leiomyo_sign <- deg_ups_leiomyo_sign[deg_ups_leiomyo_sign$logFC>0.0,]
go_enrich2 <- enrichGO(gene =   deg_ups_leiomyo_sign$...1,
                      universe = deg_ups_leiomyo$...1,
                      OrgDb = org.Hs.eg.db,
                      keyType="SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.05,
                      readable = TRUE, 
                      minGSSize = 3)

barplot(go_enrich, title = "Over Represented Pathways UPS vs LMS",showCategory = 20)
cnetplot(go_enrich,max.overlaps =200)



deg_leiomyo_lipo_sign <- deg_leiomyo_lipo[deg_leiomyo_lipo$adj.P.Val<0.05,]
deg_leiomyo_lipo_sign <- deg_leiomyo_lipo_sign[deg_leiomyo_lipo_sign$logFC>0.0,]
go_enrich3 <- enrichGO(gene =   deg_leiomyo_lipo_sign$...1,
                      universe = deg_leiomyo_lipo$...1,
                      OrgDb = org.Hs.eg.db,
                      keyType="SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 0.05,
                      readable = TRUE, 
                      minGSSize = 3)

barplot(go_enrich, title = "Over Represented Pathways LMS vs DDLPS",showCategory = 20)
cnetplot(go_enrich,max.overlaps =200)


ups_lipo <- go_enrich1@result$Description
ups_leiomyo <- go_enrich2@result$Description
lms_ddlps <- go_enrich3@result$Description

# Load library
library(VennDiagram)



# Extract GO term descriptions
ups_lipo <- data.frame(go_enrich1@result)
ups_lipo <- ups_lipo[ups_lipo$p.adjust <0.05,]$Description
ups_leiomyo <- data.frame(go_enrich2@result)
ups_leiomyo <- ups_leiomyo[ups_leiomyo$p.adjust <0.05,]$Description
lms_ddlps <- data.frame(go_enrich3@result)
lms_ddlps <- lms_ddlps[lms_ddlps$p.adjust <0.05,]$Description


# Create named list for Venn
go_lists <- list(
  "UPS vs DDLPS" = ups_lipo,
  "UPS vs LMS" = ups_leiomyo,
  "LMS vs DDLPS" = lms_ddlps
)

venn.plot <- venn.diagram(
  x = go_lists,
  filename = NULL,
  category.names = c("UPS vs DDLPS", "UPS vs LMS", "LMS vs DDLPS"),
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.4,
  cat.cex = 1.2,
  cat.pos = c(-20, 20, 180),  # adjust positions of category names
  cat.dist = c(0.05, 0.05, 0.05),
  main = "Overlapping Enriched Pathways",
  main.cex = 1.5
)

grid::grid.draw(venn.plot)

# Render the plot
grid::grid.draw(venn.plot)




