## ----setup, include=FALSE------------------------------------------------------------------------------
## directories and paths
# projectPath local
projectPath <- "/home/diego/Documentos/proyectos_master/proyectos_transcriptomica/Epigenomica/Epigenomics/Epigenomica/"
analysisPath <- file.path(projectPath, "Enriquecimiento_Funcional")
outputPath <- file.path(analysisPath, "EnrichmentSummary/annotatrGenes")


# data path github
dataPath <- "/home/diego/Documentos/proyectos_master/proyectos_transcriptomica/Transcriptomics/Epigenomica/"

dirO.7 <- "bed_E1_filtrado_0.7"
setwd(analysisPath)

if (!dir.exists(outputPath) | !dir.exists(file.path(outputPath, "Plots"))) {
  dir.create(outputPath, recursive = T)
  dir.create(file.path(outputPath, "Plots"), showWarnings = F)
}

# to files
prefix <- "EnrichmentAnalysisSum"

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
library(annotatr)
library(ggplot2)



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

barplot(ggo1$MF, showCategory = 15, font.size = 9) + ggtitle("GO Molecular Function") + 
  theme(plot.title = element_text(face = "bold"))

ChIPseeker::upsetplot(ggo1$MF)
  # upsetplot(ggo1$MF)

# goplot(ggo1)

## ------------------------------------------------------------------------------------------------------
barplot(ggo1$CC, showCategory = 10, font.size = 9) + ggtitle("GO Cellular Component") + 
  theme(plot.title = element_text(face = "bold"))


## ---- fig.asp=0.6, fig.width=10------------------------------------------------------------------------
barplot(ggo1$BP, showCategory = 20, font.size = 9) + ggtitle("GO Biological Process") + 
  theme(plot.title = element_text(face = "bold"))


## ----simplifyGO----------------------------------------------------------------------------------------
ggo1MFsimplified <- simplify(ggo1$MF)

dotplot(ggo1MFsimplified, showCategory = 20, font.size = 8)

barplot(ggo1MFsimplified, showCategory = 20, font.size = 8) + ggtitle("GO Molecular Function") + 
  theme(plot.title = element_text(face = "bold"))



## ------------------------------------------------------------------------------------------------------
kegg1 <- enrichKEGG(gene = allGenesDF$ENTREZID,
                    organism = "hsa")


## ------------------------------------------------------------------------------------------------------
dotplot(kegg1, showCategory = 20, font.size = 9) + ggtitle("KEGG pathways") + 
  theme(plot.title = element_text(face = "bold"))


barplot(kegg1, showCategory = 20, font.size = 9) + ggtitle("KEGG pathways") + 
  theme(plot.title = element_text(face = "bold"))

ChIPseeker::upsetplot(kegg1) + ggtitle("KEGG pathways") + 
  theme(plot.title = element_text(face = "bold"))


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
barplot(typeCells, "GeneRatio", showCategory = 10, font.size = 9) + ggtitle("Enriquecimiento CellMarker") + 
  theme(plot.title = element_text(face = "bold"))
ChIPseeker::upsetplot(typeCells)


## ------------------------------------------------------------------------------------------------------
  motifs <- msigdbr(species = "Homo sapiens", category = "C3") 



## ------------------------------------------------------------------------------------------------------
geneImmSginature <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)

immCells <- enricher(unique(allGenesDF$ENTREZID), TERM2GENE = geneImmSginature)


## ---- fig.asp=0.4, fig.width=9-------------------------------------------------------------------------
dotplot(immCells, showCategory = 10, font.size = 6)
barplot(immCells, "GeneRatio", showCategory = 10, font.size = 8) + ggtitle("Enriquecimiento MSigDB (C7)") + 
  theme(plot.title = element_text(face = "bold"))
ChIPseeker::upsetplot(immCells)


## ------------------------------------------------------------------------------------------------------
dso1 <- enrichDGN(allGenesDF$ENTREZID, readable = TRUE)


## ------------------------------------------------------------------------------------------------------
dotplot(dso1, showCategory = 20, font.size = 8)

barplot(dso1, showCategory = 20, font.size = 9) + ggtitle("Enriquecimiento MSigDB (C7)") + 
  theme(plot.title = element_text(face = "bold"))

ChIPseeker::upsetplot(dso1)

