###############################################################################

# Programa: GO enrichment analysis for group AS and E.
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


#### Alternative Splicing and Expression related ####
AS_y_E_df <- Datos_limpios %>% 
  dplyr::filter(Pv_wilcox < 0.05,
                P.Value < 0.05)

Lista_genes_AS_y_E <- unique(AS_y_E_df$Gene)

AS_y_E_ENTREZ <- mapIds(org.Hs.eg.db, Lista_genes_AS_y_E, 'ENTREZID', 'SYMBOL')

# GO:Celular Component - GO classification
gg_CC_AS_y_E <- groupGO(gene = AS_y_E_ENTREZ,
                        OrgDb = org.Hs.eg.db,
                        ont = "CC",
                        level = 3,
                        readable = TRUE)

head(gg_CC_AS_y_E)
barplot(gg_CC_AS_y_E, drop=TRUE, showCategory=25, vertex.label.cex=0.8, font = 6)

# GO:Celular Component - Go Over-representation analysis
go_overrepresentation_CC_AS_y_E <- enrichGO(gene = AS_y_E_ENTREZ,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "CC",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.2,
                                            readable = TRUE)
head(go_overrepresentation_CC_AS_y_E)
barplot(go_overrepresentation_CC_AS_y_E)
goplot(go_overrepresentation_CC_AS_y_E)

# GO:Molecular Function - GO classification
gg_MF_AS_y_E <- groupGO(gene = AS_y_E_ENTREZ,
                        OrgDb = org.Hs.eg.db,
                        ont = "MF",
                        level = 3,
                        readable = TRUE)

head(gg_MF_AS_y_E)
barplot(gg_MF_AS_y_E, drop=TRUE, showCategory=25, vertex.label.cex=0.8, font = 6)

# GO:Molecular Function - Go Over-representation analysis
go_overrepresentation_MF_AS_y_E <- enrichGO(gene = AS_y_E_ENTREZ,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "MF",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.2,
                                            readable = TRUE)
head(go_overrepresentation_MF_AS_y_E)
barplot(go_overrepresentation_MF_AS_y_E)
goplot(go_overrepresentation_MF_AS_y_E)

# GO:Biological Process - GO classification
gg_BP_AS_y_E <- groupGO(gene = AS_y_E_ENTREZ,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        level = 3,
                        readable = TRUE)

head(gg_BP_AS_y_E)
barplot(gg_BP_AS_y_E, drop=TRUE, showCategory=25, vertex.label.cex=0.8, font = 6)

# GO:Biological Process - Go Over-representation analysis
go_overrepresentation_BP_AS_y_E <- enrichGO(gene = AS_y_E_ENTREZ,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.2,
                                            readable = TRUE)
head(go_overrepresentation_BP_AS_y_E)
barplot(go_overrepresentation_BP_AS_y_E)
goplot(go_overrepresentation_BP_AS_y_E)

############################# FIN DEL PROGRAMA ################################