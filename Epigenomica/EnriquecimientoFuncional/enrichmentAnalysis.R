## ----setup, include=FALSE------------------------------------------------------------------------------
## directories and paths
# projectPath local
projectPath <- "/home/diego/Documentos/proyectos_master/proyectos_transcriptomica/Epigenomica/Epigenomics/Epigenomica/"
analysisPath <- file.path(projectPath, "Enriquecimiento_Funcional")
outputPath <- file.path(analysisPath, "EnrichmentAnalysis")
# data path github
dataPath <- "/home/diego/Documentos/proyectos_master/proyectos_transcriptomica/Transcriptomics/Epigenomica/"

dirO.7 <- "bed_E1_filtrado_0.7"
setwd(analysisPath)

if (!dir.exists(outputPath) | !dir.exists(file.path(outputPath, "Plots"))) {
  dir.create(outputPath, recursive = T)
  dir.create(file.path(outputPath, "Plots"), showWarnings = F)
}

# prefix to files
prefix <- "EnrichmentAnalysis"


options(stringsAsFactors = FALSE)

## packages used
library(dplyr)
library(clusterProfiler)
library(GOSemSim)
library(DOSE)
library(ReactomePA)
library(msigdbr)
library(KEGGprofile)
library(org.Hs.eg.db)

library(ChIPpeakAnno)
library(annotatr)
library(ChIPseeker)

# library(TxDb.Hsapiens.UCSC.hg19.knownGene)



## ----dataLoading---------------------------------------------------------------------------------------
dataAnnotatr <- read.csv(file.path(dataPath, dirO.7, "annotatr_results_E1_0.7.csv"))


## ----allGenes------------------------------------------------------------------------------------------
allGenesDF <- data.frame(ENTREZID = dataAnnotatr$annot.gene_id,
                       SYMBOL = dataAnnotatr$annot.symbol)

# allGenesDF <- bitr(allGenes, fromType = "SYMBOL",
#                       toType = c("ENSEMBL", "ENTREZID"),
#                       OrgDb = org.Hs.eg.db)


## ----goAllGenes----------------------------------------------------------------------------------------
ggo1 <- sapply(c("MF", "CC", "BP"), function(term) {
              enrichGO(unique(allGenesDF$ENTREZID), org.Hs.eg.db,
                             ont = term, qvalueCutoff = 0.05,
                             readable = TRUE)})


## ----plotGoAllGenes------------------------------------------------------------------------------------
dotplot(ggo1$MF, showCategory = 20, font.size = 8)

barplot(ggo1$MF, showCategory = 20, font.size = 8)

upsetplot(ggo1$MF)

# goplot(ggo1)


## ----goTableAllGenes-----------------------------------------------------------------------------------
DT::datatable(ggo1$MF@result)


## ------------------------------------------------------------------------------------------------------
ggo1MFsimplified <- simplify(ggo1$MF)

dotplot(ggo1MFsimplified, showCategory = 20, font.size = 8)

barplot(ggo1MFsimplified, showCategory = 20, font.size = 8)


## ------------------------------------------------------------------------------------------------------
dataHumanCells <- read.delim(file.path(analysisPath, "/data/Human_cell_markers.txt"))

dataHumanCells <- dataHumanCells %>% tidyr::unite("cellMarker", tissueType, 
                                                  cancerType, cellName, sep=", ") %>% 
   dplyr::select(cellMarker, geneID) %>%
   dplyr::mutate(geneID = strsplit(geneID, ', '))


## ------------------------------------------------------------------------------------------------------
typeCells <- enricher(unique(allGenesDF$ENTREZID), TERM2GENE = dataHumanCells,
                      pAdjustMethod = "fdr", minGSSize = 5, qvalueCutoff = 0.2)

dotplot(typeCells, showCategory = 10, font.size = 8)
barplot(typeCells, "GeneRatio", showCategory = 10, font.size = 8)
upsetplot(typeCells)


## ------------------------------------------------------------------------------------------------------
DT::datatable(typeCells@result)


## ------------------------------------------------------------------------------------------------------
dataBloodCells <- read.delim(file.path(analysisPath, "/data/Human_cell_markers.txt"))

monocytesInfo <- unique(dataBloodCells[grepl(".*([mM]onocyte)", dataBloodCells$cellName),]$tissueType)

dataBloodCells <- dataBloodCells %>% filter(tissueType %in% monocytesInfo) %>% 
                    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
                    dplyr::select(cellMarker, geneID) %>%
                    dplyr::mutate(geneID = strsplit(geneID, ', '))


## ------------------------------------------------------------------------------------------------------
typeCellsBlood <- enricher(unique(allGenesDF$ENTREZID), TERM2GENE = dataBloodCells,
                      pAdjustMethod = "fdr", minGSSize = 5, qvalueCutoff = 0.2)

dotplot(typeCellsBlood, showCategory = 10, font.size = 8)
barplot(typeCellsBlood, "GeneRatio", showCategory = 10, font.size = 8)


## ------------------------------------------------------------------------------------------------------
DT::datatable(typeCellsBlood@result)


## ------------------------------------------------------------------------------------------------------
geneImmSginature <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)

immCells <- enricher(unique(allGenesDF$ENTREZID), TERM2GENE = geneImmSginature)


## ------------------------------------------------------------------------------------------------------
dotplot(immCells, showCategory = 10, font.size = 6)
barplot(immCells, "GeneRatio", showCategory = 10, font.size = 6)
upsetplot(immCells)


## ------------------------------------------------------------------------------------------------------
dso1 <- enrichDGN(allGenesDF$ENTREZID, readable = TRUE)


## ------------------------------------------------------------------------------------------------------
dotplot(dso1, showCategory = 20, font.size = 8)

barplot(dso1, showCategory = 20, font.size = 8)

upsetplot(dso1)

## ------------------------------------------------------------------------------------------------------
do1 <- enrichDO(allGenesDF$ENTREZID,
                 readable = TRUE)


## ------------------------------------------------------------------------------------------------------
dotplot(do1, showCategory = 20, font.size = 8)

barplot(do1, showCategory = 20, font.size = 8)

upsetplot(do1)


## ------------------------------------------------------------------------------------------------------
kegg1 <- enrichKEGG(gene = allGenesDF$ENTREZID,
                    organism = "hsa")


## ------------------------------------------------------------------------------------------------------
dotplot(kegg1, showCategory = 20, font.size = 8)

barplot(kegg1, showCategory = 20, font.size = 8)

upsetplot(kegg1)


## ------------------------------------------------------------------------------------------------------
pho_KEGGresult <- find_enriched_pathway(unique(allGenesDF$ENTREZID), 
                                        species = 'hsa')

DT::datatable(pho_KEGGresult[[1]])


## ----promoterGenes-------------------------------------------------------------------------------------
genesPromoter <- dataAnnotatr %>% filter(grepl("^(promoter).*", annot.id),
                                         !is.na(annot.gene_id))
genesPromoterDF <- data.frame(ENTREZID = genesPromoter$annot.gene_id,
                              SYMBOL = genesPromoter$annot.symbol)
 
# genesPromoterDF <- bitr(genesPromoter, fromType = "SYMBOL",
#                       toType = c("ENSEMBL", "ENTREZID"),
#                       OrgDb = org.Hs.eg.db)


## ----goPromoterGenes-----------------------------------------------------------------------------------
ggo2 <- sapply(c(c("MF", "CC", "BP")), function(x) {
  enrichGO(unique(genesPromoterDF$ENTREZID), org.Hs.eg.db,
                 ont = "MF", qvalueCutoff = 0.05,
                 readable = TRUE)}
  )


## ----plotGoPromoterGenes-------------------------------------------------------------------------------
dotplot(ggo2$MF, showCategory = 20, font.size = 8)

barplot(ggo2$MF, showCategory = 20, font.size = 8)

upsetplot(ggo2$MF)

goplot(ggo2$MF)


## ----tablePromoterGenes--------------------------------------------------------------------------------
DT::datatable(ggo2$MF@result)


## ------------------------------------------------------------------------------------------------------
ggo2MFsimplified <- simplify(ggo2$MF)

dotplot(ggo2MFsimplified, showCategory = 20, font.size = 8)

barplot(ggo2MFsimplified, showCategory = 20, font.size = 8)


## ------------------------------------------------------------------------------------------------------
hsGO <- godata('org.Hs.eg.db', ont="MF")

distan <- mgeneSim(head(genesPromoterDF$ENTREZID, n = 20), semData = hsGO, 
                   measure = "Wang", combine = "BMA")

pheatmap::pheatmap(distan)


## ------------------------------------------------------------------------------------------------------
typeCellsPromoters <- enricher(unique(genesPromoterDF$ENTREZID), TERM2GENE = dataHumanCells,
                      pAdjustMethod = "fdr", minGSSize = 5, qvalueCutoff = 0.2)

dotplot(typeCellsPromoters, showCategory = 50, font.size = 8)
barplot(typeCellsPromoters, showCategory = 50, font.size = 8)


## ------------------------------------------------------------------------------------------------------
DT::datatable(typeCellsPromoters@result)


## ------------------------------------------------------------------------------------------------------
dso2 <- enrichDGN(genesPromoterDF$ENTREZID,
                 readable = TRUE)


## ------------------------------------------------------------------------------------------------------
dotplot(dso2, showCategory = 20, font.size = 8)

barplot(dso2, showCategory = 20, font.size = 8)

upsetplot(dso2)


## ------------------------------------------------------------------------------------------------------
enrichplot::pmcplot(head(dso2$Description), 2012:2019)


## ----annotationBED, warning=FALSE, message=FALSE-------------------------------------------------------

if (file.exists(file.path(outputPath, "dm_annotated.rds")) &
    file.exists(file.path(outputPath, "dm_random_annotated.rds"))){
  
  dm_annotated <- readRDS(file.path(outputPath, "dm_annotated.rds"))
  dm_random_annotated <- readRDS(file.path(outputPath, "dm_random_annotated.rds"))
  df_dm_annotated <- data.frame(dm_annotated)
  
} else {
  
  regionsInter <- read_regions(con = file.path(analysisPath, "data/intersect_iPScells/solo_monocitos.bed"), 
                             genome = 'hg19',
                             format = 'bed')

  annots <- c('hg19_genes_cds', 'hg19_basicgenes', 'hg19_cpgs', 'hg19_genes_intergenic',
             'hg19_enhancers_fantom')
  
  annotations <- build_annotations(genome = 'hg19', annotations = annots)
  # Intersect the regions we read in with the annotations
  dm_annotated <- annotate_regions(regions = regionsInter, 
                                   annotations = annotations,
                                   minoverlap = 100L, #Overlap of annotations, f=-0.5
                                   ignore.strand = TRUE,
                                   quiet = FALSE)
  
  df_dm_annotated <- data.frame(dm_annotated)
  
  # Randomize the input regions
  dm_random_regions <- randomize_regions(regions = regionsInter,
                                        allow.overlaps = TRUE,
                                        per.chromosome = TRUE)
  
  # Annotate the random regions using the same annotations as above
  # These will be used in later functions
  dm_random_annotated <- annotate_regions(regions = dm_random_regions,
                                          annotations = annotations,
                                          minoverlap = 100L, #Overlap of annotations, f=-0.5
                                          ignore.strand = TRUE, quiet = TRUE)
  
  saveRDS(dm_annotated, file.path(outputPath, "dm_annotated.rds"))
  saveRDS(dm_random_annotated, file.path(outputPath, "dm_random_annotated.rds"))
}




## ------------------------------------------------------------------------------------------------------
# Find the number of regions per annotation type
dm_annsum <- summarize_annotations(annotated_regions = dm_annotated, quiet = TRUE)
print(dm_annsum)

# Count the occurrences of classifications in the Status
# column across the annotation types.
dm_catsum <- summarize_categorical(annotated_regions = dm_annotated,
                                  quiet = TRUE)
print(dm_catsum)


## ------------------------------------------------------------------------------------------------------
#View a heatmap of regions occurring in pairs of annotations
annots_order <- c('hg19_genes_promoters', 'hg19_genes_5UTRs', 
                  'hg19_genes_exons', 'hg19_genes_introns',
                  'hg19_genes_3UTRs', 'hg19_genes_cds',
                  'hg19_enhancers_fantom')

dm_vs_coannotations <- plot_coannotations(annotated_regions = dm_annotated,
                                          annotation_order = annots_order, 
                                          axes_label = 'Annotations',
                                          plot_title = 'Regions in Pairs of Annotations')
print(dm_vs_coannotations)

dm_annotations_plot <- plot_annotation(annotated_regions = dm_annotated,
                annotation_order = annots_order,
                plot_title = 'Annotations for monocytes-iPCs',
                x_label = 'Annotation type',
                y_label = 'Count')

print(dm_annotations_plot)


## ------------------------------------------------------------------------------------------------------
# View the number of regions per annotation and include the annotation
# of randomized regions
dm_annotations_plot_wrandom = plot_annotation(annotated_regions = dm_annotated,
                                              annotated_random = dm_random_annotated,
                                              annotation_order = annots_order,
                                              plot_title = 'Annotations for monocytes-iPC (with rndm.)',
                                              x_label = 'Annotation type',
                                              y_label = 'Count')
print(dm_annotations_plot_wrandom)

#CpGIslands
annots_order = c('hg19_cpg_islands', 'hg19_cpg_shores',
                 'hg19_cpg_shelves', 'hg19_cpg_inter')

dm_annotations_plot_wrandom = plot_annotation(annotated_regions = dm_annotated,
                                              annotated_random = dm_random_annotated,
                                              annotation_order = annots_order,
                                              plot_title = 'Annotations for monocytes-iPC (with rndm.)', 
                                              x_label = 'Annotation type',
                                              y_label = 'Count')
print(dm_annotations_plot_wrandom)




## ------------------------------------------------------------------------------------------------------
genesMonoIPC <- data.frame(ENTREZID = df_dm_annotated$annot.gene_id,
                           SYMBOL = df_dm_annotated$annot.symbol)

intersectGenes <- intersect(genesMonoIPC$ENTREZID, allGenesDF$ENTREZID)

cat(">> Número de genes en E1:", length(unique(allGenesDF$ENTREZID)))
cat("\n>> Número de genes en Mono-iPCs:", length(unique(genesMonoIPC$ENTREZID)))
cat("\n>> Número de genes coincidentes:", length(intersectGenes))


## ------------------------------------------------------------------------------------------------------
ggo3 <- sapply(c("MF", "CC", "BP"), function(term) {
              enrichGO(unique(genesMonoIPC$ENTREZID), org.Hs.eg.db,
                             ont = term, qvalueCutoff = 0.05,
                             readable = TRUE)})


## ------------------------------------------------------------------------------------------------------
dotplot(ggo3$MF, showCategory = 20, font.size = 8)

barplot(ggo3$MF, showCategory = 20, font.size = 8)

upsetplot(ggo3$MF)

# goplot(ggo1)



## ------------------------------------------------------------------------------------------------------
kegg3 <- enrichKEGG(unique(genesMonoIPC$ENTREZID),
                    organism = "hsa")


## ------------------------------------------------------------------------------------------------------
dotplot(kegg3, showCategory = 20, font.size = 8)

barplot(kegg3, showCategory = 20, font.size = 8)

upsetplot(kegg3)

