library(clusterProfiler)
library(ggplot2)
library(plotly)
library(msigdbr)
library(dplyr)
library(gridExtra)
library(org.Hs.eg.db)
library(rempsyc)

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

hs_kegg_df <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat == "CP:KEGG" # This is because we only want KEGG pathways
  ) 
hs_kegg_df_with_gs_source <- hs_kegg_df[,c("gs_name","gene_symbol","gs_exact_source")] %>%subset(grepl("^hsa00|^hsa01|^hsa02", gs_exact_source))


hs_kegg_df <- hs_kegg_df[,c("gs_name","gene_symbol")]
hs_kegg_df_1 <- hs_kegg_df
# filtrar genes que aparecem repetidos
hs_kegg_df <- hs_kegg_df %>%
  distinct(gs_name, gene_symbol, .keep_all = TRUE)

# filtrar genes repetidos em cada pathway!
gene_names <-row.names(dge_leiomyo_lipo)
gene_names_2 <- row.names(dge_ups_leiomyo)
gene_names_3 <- row.names(dge_ups_lipo)
genes <- row.names(dge_leiomyo_lipo_over)
genes2 <- row.names(dge_ups_leiomyo_over)
genes3 <- row.names(dge_ups_lipo_over)

genes_under <- row.names(dge_leiomyo_lipo_under)
genes2_under <- row.names(dge_ups_leiomyo_under)
genes3_under <- row.names(dge_ups_lipo_under)

#gene symbol to all get pathways with genes in common
matching_genes <- hs_kegg_df$gs_name[hs_kegg_df$gene_symbol %in% gene_names]

keep_names <- unique(matching_genes)

pathway_gene_counts <- hs_kegg_df %>%
  group_by(gs_name) %>%
  summarise(number_of_genes_involved = n_distinct(gene_symbol))

# #keep only the pathways that have at least one matching gene with kit
hs_kegg_df_kit <- hs_kegg_df[hs_kegg_df$gs_name %in% keep_names,]
hs_kegg_df <- hs_kegg_df[hs_kegg_df$gs_name %in% keep_names,]
hs_kegg_df_with_gs_source <- hs_kegg_df_with_gs_source[hs_kegg_df_with_gs_source$gs_name %in% keep_names,]

#keep only genes involved in kit
hs_kegg_df_kit <- hs_kegg_df_kit[hs_kegg_df_kit$gene_symbol %in% gene_names,]
hs_kegg_df_with_gs_source <- hs_kegg_df_with_gs_source[hs_kegg_df_with_gs_source$gene_symbol %in% gene_names,]
hs_kegg_df_kit <-hs_kegg_df_kit[hs_kegg_df_kit$gs_name %in% hs_kegg_df_with_gs_source$gs_name,]

write.csv(hs_kegg_df_kit,"RESULTS/METABOLIC_GENES_IN_KIT.csv")


pathway_gene_counts_included_in_kit <- hs_kegg_df_kit %>%
  group_by(gs_name) %>%
  summarise(number_of_genes_involved = n_distinct(gene_symbol))


#genes dge leiomyo lipo
#DGE
hs_kegg_df_kit_1 <- hs_kegg_df_kit[hs_kegg_df_kit$gene_symbol %in% genes,]
pathway_gene_counts_included_in_kit_1 <- hs_kegg_df_kit_1 %>%
  group_by(gs_name) %>%
  summarise(number_of_genes_involved = n_distinct(gene_symbol))
#UNDER
hs_kegg_df_kit_1_under <- hs_kegg_df_kit[hs_kegg_df_kit$gene_symbol %in% genes_under,]
pathway_gene_counts_included_in_kit_1_under <- hs_kegg_df_kit_1_under %>%
  group_by(gs_name) %>%
  summarise(number_of_genes_involved = n_distinct(gene_symbol))


#genes dge nos leiomyo
hs_kegg_df_kit_2 <- hs_kegg_df_kit[hs_kegg_df_kit$gene_symbol %in% genes2,]
pathway_gene_counts_included_in_kit_2 <- hs_kegg_df_kit_2 %>%
  group_by(gs_name) %>%
  summarise(number_of_genes_involved = n_distinct(gene_symbol))

#UNDER
hs_kegg_df_kit_2_under <- hs_kegg_df_kit[hs_kegg_df_kit$gene_symbol %in% genes2_under,]
pathway_gene_counts_included_in_kit_2_under <- hs_kegg_df_kit_2_under %>%
  group_by(gs_name) %>%
  summarise(number_of_genes_involved = n_distinct(gene_symbol))

#genes dge nos lipo
hs_kegg_df_kit_3 <- hs_kegg_df_kit[hs_kegg_df_kit$gene_symbol %in% genes3,]
pathway_gene_counts_included_in_kit_3 <- hs_kegg_df_kit_3 %>%
  group_by(gs_name) %>%
  summarise(number_of_genes_involved = n_distinct(gene_symbol))

#UNDER
hs_kegg_df_kit_3_under <- hs_kegg_df_kit[hs_kegg_df_kit$gene_symbol %in% genes3_under,]
pathway_gene_counts_included_in_kit_3_under <- hs_kegg_df_kit_3_under %>%
  group_by(gs_name) %>%
  summarise(number_of_genes_involved = n_distinct(gene_symbol))


hs_kegg_df_kit_1 <- hs_kegg_df_kit_1 %>%
  group_by(gs_name) %>%
  summarize(genes_leiomyo_lipo = paste(gene_symbol, collapse = "-"))

hs_kegg_df_kit_2 <- hs_kegg_df_kit_2 %>%
  group_by(gs_name) %>%
  summarize(genes_nos_leiomyo = paste(gene_symbol, collapse = "-"))

hs_kegg_df_kit_3 <- hs_kegg_df_kit_3 %>%
  group_by(gs_name) %>%
  summarize(genes_nos_lipo = paste(gene_symbol, collapse = "-"))

#under
hs_kegg_df_kit_1_under <- hs_kegg_df_kit_1_under %>%
  group_by(gs_name) %>%
  summarize(genes_leiomyo_lipo = paste(gene_symbol, collapse = "-"))

hs_kegg_df_kit_2_under <- hs_kegg_df_kit_2_under %>%
  group_by(gs_name) %>%
  summarize(genes_nos_leiomyo = paste(gene_symbol, collapse = "-"))

hs_kegg_df_kit_3_under <- hs_kegg_df_kit_3_under %>%
  group_by(gs_name) %>%
  summarize(genes_nos_lipo = paste(gene_symbol, collapse = "-"))

#
final_result <- left_join(pathway_gene_counts, pathway_gene_counts_included_in_kit, by = "gs_name")

final_result <- left_join(final_result, pathway_gene_counts_included_in_kit_1,  by = "gs_name")

final_result <- left_join(final_result, pathway_gene_counts_included_in_kit_2,  by = "gs_name")

final_result <- left_join(final_result, pathway_gene_counts_included_in_kit_3,  by = "gs_name")

final_result <- left_join(final_result, hs_kegg_df_kit_1,  by = "gs_name")

final_result <- left_join(final_result, hs_kegg_df_kit_2,  by = "gs_name")
final_result <- left_join(final_result, hs_kegg_df_kit_3,  by = "gs_name")

final_result <- left_join(final_result, hs_kegg_df_kit_1_under,  by = "gs_name")

final_result <- left_join(final_result, hs_kegg_df_kit_2_under,  by = "gs_name")

final_result <- left_join(final_result, hs_kegg_df_kit_3_under,  by = "gs_name")

colnames(final_result) <- c("pathway", "genes_involved_in_pathway", "kit_gene_counts",
                            "leiomyo_lipo_gene_counts","nos_leiomyo_gene_counts", "nos_lipo_gene_counts",
                            "genes_leiomyo_lipo","genes_nos_leiomyo","genes_nos_lipo","genes_leiomyo_lipo_commonly_expressed"
                            ,"genes_nos_leiomyo_commonly_expressed","genes_nos_lipo_commonly_expressed")

# Replace NA with 0 in the kit_gene_counts column
final_result$kit_gene_counts[is.na(final_result$kit_gene_counts)] <- 0

final_result$leiomyo_lipo_gene_counts[is.na(final_result$leiomyo_lipo_gene_counts)] <- 0

final_result$nos_leiomyo_gene_counts[is.na(final_result$nos_leiomyo_gene_counts)] <- 0

final_result$nos_lipo_gene_counts[is.na(final_result$nos_lipo_gene_counts)] <- 0

final_result <- final_result[final_result$kit_gene_counts != 0,]

final_result <- final_result[order(final_result$genes_involved_in_pathway, decreasing = TRUE),]


final_result <- final_result[final_result$kit_gene_counts > 1 ,]

ggplot(final_result, aes(x = reorder(pathway, -genes_involved_in_pathway))) +
  geom_bar(aes(y = genes_involved_in_pathway), stat = "identity", fill = "blue", width = 0.5) +
  geom_bar(aes(y = kit_gene_counts), stat = "identity", fill = "orange", width = 0.5) +
  geom_text(aes(y = genes_involved_in_pathway, label = genes_involved_in_pathway), vjust = -0.2, color = "black") +
  geom_text(aes(y = kit_gene_counts, label = kit_gene_counts), vjust = -0.2, color = "white") +
  labs(title = "Concordance between Genes in Pathways and Genes in Kit", x = "Pathway", y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#
colnames(final_result)


gg <- ggplot(final_result, aes(x = reorder(pathway, -kit_gene_counts), text = paste("Normal:", genes_leiomyo_lipo_commonly_expressed," DGE: ",genes_leiomyo_lipo))) +
  geom_bar(aes(y = kit_gene_counts), stat = "identity", fill = "lightblue", width = 0.5) +
  geom_bar(aes(y = leiomyo_lipo_gene_counts), stat = "identity", fill = "orange", width = 0.5) +
  geom_text(aes(y = kit_gene_counts, label = kit_gene_counts), vjust = -0.2, color = "black") +
  geom_text(aes(y = leiomyo_lipo_gene_counts, label = leiomyo_lipo_gene_counts), vjust = -0.2, color = "black") +
  labs(title = "Concordance between Genes from Kit in Pathways and DGE genes Leiomyo Lipo", x = "Pathway", y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

ggplotly(gg)

gg2 <- ggplot(final_result, aes(x = reorder(pathway, -kit_gene_counts), text = paste("Normal:", genes_nos_leiomyo_commonly_expressed," DGE: ",genes_nos_leiomyo))) +
  geom_bar(aes(y = kit_gene_counts), stat = "identity", fill = "lightblue", width = 0.5) +
  geom_bar(aes(y = nos_leiomyo_gene_counts), stat = "identity", fill = "orange", width = 0.5) +
  geom_text(aes(y = kit_gene_counts, label = kit_gene_counts), vjust = -0.2, color = "black") +
  geom_text(aes(y = nos_leiomyo_gene_counts, label = nos_leiomyo_gene_counts), vjust = -0.2, color = "black") +
  labs(title = "Concordance between Genes from Kit in Pathways and DGE genes NOS Leiomyo", x = "Pathway", y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplotly(gg2)

gg3 <- ggplot(final_result, aes(x = reorder(pathway, -kit_gene_counts), text = paste("Normal:", genes_nos_lipo_commonly_expressed," DGE: ",genes_nos_lipo))) +
  geom_bar(aes(y = kit_gene_counts), stat = "identity", fill = "lightblue", width = 0.5) +
  geom_bar(aes(y = nos_lipo_gene_counts), stat = "identity", fill = "orange", width = 0.5) +
  geom_text(aes(y = kit_gene_counts, label = kit_gene_counts), vjust = -0.2, color = "black") +
  geom_text(aes(y = nos_lipo_gene_counts, label = nos_lipo_gene_counts), vjust = -0.2, color = "black") +
  labs(title = "Concordance between Genes from Kit in Pathways and DGE genes NOS Lipo", x = "Pathway", y = "Number of Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplotly(gg3)


grid.arrange(gg,gg2,gg3,ncol =3)

write.csv(final_result, "RESULTS/Metabolic_Pathway_information.csv")
nice_table(
  final_result
)

go_bp_genesets = msigdbr(species = "human", category = "C5", subcategory = "BP")
go_bp_genesets = go_bp_genesets[go_bp_genesets$gene_symbol %in% hs_kegg_df_kit$gene_symbol,]

go_bp_genesets = go_bp_genesets[go_bp_genesets$gene_symbol %in% c("SDHA","SDHB","SDHC","SDHD"),]

write.csv(go_bp_genesets,"RESULTS/biological_processes_with_SDH.csv")
go_mf_genesets = msigdbr(species = "human", category = "C5", subcategory = "MF")
go_mf_genesets = go_mf_genesets[go_mf_genesets$gene_symbol %in% hs_kegg_df_kit$gene_symbol,]
go_mf_genesets = go_mf_genesets[go_mf_genesets$gene_symbol %in% c("SDHA","SDHB","SDHC","SDHD"),]
write.csv(go_mf_genesets,"RESULTS/molecular_functions_with_SDH.csv")

go_cc_genesets = msigdbr(species = "human", category = "C5", subcategory = "CC")
go_cc_genesets = go_cc_genesets[go_cc_genesets$gene_symbol %in% hs_kegg_df_kit$gene_symbol,]
go_cc_genesets = go_cc_genesets[go_cc_genesets$gene_symbol %in% c("SDHA","SDHB","SDHC","SDHD"),]
write.csv(go_cc_genesets,"RESULTS/celullar_components_with_SDH.csv")



