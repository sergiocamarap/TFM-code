###############################################################################

# Programa: GO enrichment analysis for group AS related.
# Hecho por: Sergio Cámara Peña
# Fecha: 23/06/2022
# Versión: 1.0
# Encoding: UTF-8

# Nota: Para ver más información acerca de los grupos ir al R Markdown.

###############################################################################

#### Cargamos los paquetes necesarios para el programa ####
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)


#### Lectura y limpiza de datos ####
Datos_originales <- read_delim("Documentos/20222_01_21_Hit_wilcox_expr_Surv.csv", 
                               delim = ";", escape_double = FALSE, 
                               locale = locale(decimal_mark = ",", 
                                               grouping_mark = "."), 
                               trim_ws = TRUE)

Datos_limpios <- Datos_originales %>% 
  dplyr::select(-c("...1","Pathway","Reference"))

Lista_total_genes <- unique(Datos_limpios$Gene)


#### Related with Alternative Splicing (no focus on Expression) ####
AS_related_df <- Datos_limpios %>% 
  dplyr::filter(Pv_wilcox < 0.05)

Lista_genes_AS_related <- unique(AS_related_df$Gene)

AS_related_ENTREZ <- mapIds(org.Hs.eg.db, Lista_genes_AS_related, 'ENTREZID', 'SYMBOL')


# GO:Celular Component - GO classification
gg_CC_AS_related <- groupGO(gene = AS_related_ENTREZ,
                           OrgDb = org.Hs.eg.db,
                           ont = "CC",
                           level = 3,
                           readable = TRUE)

head(gg_CC_AS_related)
barplot(gg_CC_AS_related, drop=TRUE, showCategory=25, vertex.label.cex=0.8, font = 6)

# GO:Celular Component - Go Over-representation analysis
go_overrepresentation_CC_AS_related <- enrichGO(gene = AS_related_ENTREZ,
                                               OrgDb = org.Hs.eg.db,
                                               ont = "CC",
                                               pAdjustMethod = "BH",
                                               pvalueCutoff = 0.05,
                                               qvalueCutoff = 0.2,
                                               readable = TRUE)
head(go_overrepresentation_CC_AS_related)
barplot(go_overrepresentation_CC_AS_related)
goplot(go_overrepresentation_CC_AS_related)

# GO:Molecular Function - GO classification
gg_MF_AS_related <- groupGO(gene = AS_related_ENTREZ,
                           OrgDb = org.Hs.eg.db,
                           ont = "MF",
                           level = 3,
                           readable = TRUE)

head(gg_MF_AS_related)
barplot(gg_MF_AS_related, drop=TRUE, showCategory=25, vertex.label.cex=0.8, font = 6)

# GO:Molecular Function - Go Over-representation analysis
go_overrepresentation_MF_AS_related <- enrichGO(gene = AS_related_ENTREZ,
                                               OrgDb = org.Hs.eg.db,
                                               ont = "MF",
                                               pAdjustMethod = "BH",
                                               pvalueCutoff = 0.05,
                                               qvalueCutoff = 0.2,
                                               readable = TRUE)
head(go_overrepresentation_MF_AS_related)
barplot(go_overrepresentation_MF_AS_related)
goplot(go_overrepresentation_MF_AS_related)

# GO:Biological Process - GO classification
gg_BP_AS_related <- groupGO(gene = AS_related_ENTREZ,
                           OrgDb = org.Hs.eg.db,
                           ont = "BP",
                           level = 3,
                           readable = TRUE)

head(gg_BP_AS_related)
barplot(gg_BP_AS_related, drop=TRUE, showCategory=25, vertex.label.cex=0.8, font = 6)

# GO:Biological Process - Go Over-representation analysis
go_overrepresentation_BP_AS_related <- enrichGO(gene = AS_related_ENTREZ,
                                               OrgDb = org.Hs.eg.db,
                                               ont = "BP",
                                               pAdjustMethod = "BH",
                                               pvalueCutoff = 0.05,
                                               qvalueCutoff = 0.2,
                                               readable = TRUE)
head(go_overrepresentation_BP_AS_related)
barplot(go_overrepresentation_BP_AS_related)
goplot(go_overrepresentation_BP_AS_related)

############################# FIN DEL PROGRAMA ################################