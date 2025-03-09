library(clusterProfiler)
library(ggplot2)
library(plotly)
library(msigdbr)
library(dplyr)
library(gridExtra)
library(org.Hs.eg.db)
library(rempsyc)

dge_ups_lipo <- read.csv("RESULTS/deg_tcga_ups_lipo.csv", row.names=1)
dge_ups_lipo$genes <- NULL
names <- row.names(dge_ups_lipo)
dge_ups_lipo <- data.frame(sapply(dge_ups_lipo, as.numeric))
row.names(dge_ups_lipo) <- names


dge_ups_lipo_sign <- dge_ups_lipo[dge_ups_lipo$adj.P.Val <= 0.05,]
dge_ups_lipo_over <- dge_ups_lipo_sign[dge_ups_lipo_sign$logFC > 0,]
dge_ups_lipo_under <- dge_ups_lipo_sign[dge_ups_lipo_sign$logFC< 0,]

dge_ups_leiomyo <- read.csv("RESULTS/deg_tcga_ups_leiomyo.csv", row.names = 1)
dge_ups_leiomyo$genes <- NULL
names <- row.names(dge_ups_leiomyo)
dge_ups_leiomyo <- data.frame(sapply(dge_ups_leiomyo, as.numeric))
row.names(dge_ups_leiomyo) <- names

dge_ups_leiomyo_sign <- dge_ups_leiomyo[dge_ups_leiomyo$adj.P.Val <= 0.05,]
dge_ups_leiomyo_over <- dge_ups_leiomyo_sign[dge_ups_leiomyo_sign$logFC > 0,]
dge_ups_leiomyo_under <- dge_ups_leiomyo_sign[dge_ups_leiomyo_sign$logFC < 0,]

dge_leiomyo_lipo <- read.csv("RESULTS/deg_tcga_leiomyo_lipo.csv", row.names = 1)
dge_leiomyo_lipo$genes <- NULL
names <- row.names(dge_leiomyo_lipo)
dge_leiomyo_lipo <- data.frame(sapply(dge_leiomyo_lipo, as.numeric))
row.names(dge_leiomyo_lipo) <- names
dge_leiomyo_lipo_sign <- dge_leiomyo_lipo[dge_leiomyo_lipo$adj.P.Val <= 0.05,]
dge_leiomyo_lipo_over <- dge_leiomyo_lipo_sign[dge_leiomyo_lipo_sign$logFC > 0,]
dge_leiomyo_lipo_under <- dge_leiomyo_lipo_sign[dge_leiomyo_lipo_sign$logFC < 0,]


hs_kegg_df <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat == "CP:KEGG" # This is because we only want KEGG pathways
  )

gene_list <- dge_ups_leiomyo$logFC
names(gene_list) <- row.names(dge_ups_leiomyo)
gene_list = sort(gene_list, decreasing = TRUE)
hs_kegg_df <- hs_kegg_df[,c("gs_name","gene_symbol")]
deg2_gsea <- GSEA(geneList = gene_list,TERM2GENE = hs_kegg_df)

# Define the pathways to remove
pathways_to_remove <- c(
  "KEGG_AUTOIMUNE_THYROID_DISEASE",
  "KEGG_TYPE_I_DIABETES_MELLITUS",
  "KEGG_ASTHMA",
  "KEGG_VIRAL_MYOCARDITIS",
  "KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS",
  "KEGG_ALDOSTERONE_REGULATED_SODIUM_REABSORPTION",
  "KEGG_GLIOMA",
  "KEGG_MELANOMA",
  "KEGG_COMPLEMENT_AND_COAGULATION",
  "KEGG_AXON_GUIDANCE",
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "KEGG_PARKINSONS_DISEASE",
  "KEGG_CHRONIC_MYELOID_LEUKEMIA"
)

print(row.names(rownames(deg2_gsea@result)))

# Remove the specified pathways
data <- deg2_gsea@result[!rownames(deg2_gsea@result) %in% pathways_to_remove, ]

deg2_gsea@result <- data

ridgeplot(deg2_gsea, label_format = 50,showCategory = 15)


hs_kegg_df <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat == "CP:KEGG" # This is because we only want KEGG pathways
  )

gene_list <- dge_ups_lipo$logFC
names(gene_list) <- row.names(dge_ups_lipo)
gene_list = sort(gene_list, decreasing = TRUE)
hs_kegg_df <- hs_kegg_df[,c("gs_name","gene_symbol")]
deg2_gsea <- GSEA(geneList = gene_list,TERM2GENE = hs_kegg_df)

# Define the pathways to remove
pathways_to_remove <- c(
  "KEGG_AUTOIMUNE_THYROID_DISEASE",
  "KEGG_TYPE_I_DIABETES_MELLITUS",
  "KEGG_ASTHMA",
  "KEGG_VIRAL_MYOCARDITIS",
  "KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS",
  "KEGG_ALDOSTERONE_REGULATED_SODIUM_REABSORPTION",
  "KEGG_GLIOMA",
  "KEGG_MELANOMA",
  "KEGG_COMPLEMENT_AND_COAGULATION",
  "KEGG_AXON_GUIDANCE",
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
  "KEGG_PARKINSONS_DISEASE",
  "KEGG_CHRONIC_MYELOID_LEUKEMIA"
)

print(row.names(rownames(deg2_gsea@result)))

# Remove the specified pathways
data <- deg2_gsea@result[!rownames(deg2_gsea@result) %in% pathways_to_remove, ]


deg2_gsea@result <- data

ridgeplot(deg2_gsea, label_format = 50,showCategory = 20)


