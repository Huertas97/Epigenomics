library('annotatr')
title_str <- 'Status'
score <- 'Post.P'
#Load BED file
#Fichero de prueba solo con 2000 lineas del fichero de segmentos del Monocito 1.
dm_regions = read_regions(con = 'intersect_E1_f0.7.bed', genome = 'hg19', format = 'bed',
                          rename_name = title_str, rename_score = score)

# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots = c('hg19_genes_cds', 'hg19_basicgenes', 'hg19_cpgs', 'hg19_genes_intergenic',
           'hg19_enhancers_fantom') # cdsm genes, islas cpgs, regiones intergenicas y enhancers
# Build the annotations (a single GRanges object)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
annotations = build_annotations(genome = 'hg19', annotations = annots)
# Intersect the regions we read in with the annotations
dm_annotated = annotate_regions(
  regions = dm_regions,
  annotations = annotations,
  minoverlap = 100L, #Overlap of annotations, f=-0.5
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
print(dm_annotated)

# Coerce to a data.frame
df_dm_annotated = data.frame(dm_annotated)
# See the GRanges column of dm_annotaed expanded
print(head(df_dm_annotated))

# Randomize the input regions
dm_random_regions = randomize_regions(
  regions = dm_regions,
  allow.overlaps = TRUE,
  per.chromosome = TRUE)
# Annotate the random regions using the same annotations as above
# These will be used in later functions
dm_random_annotated = annotate_regions(
  regions = dm_random_regions,
  annotations = annotations,
  minoverlap = 100L, #Overlap of annotations, f=-0.5
  ignore.strand = TRUE,
  quiet = TRUE)

# Find the number of regions per annotation type
dm_annsum = summarize_annotations(
  annotated_regions = dm_annotated,
  quiet = TRUE)
print(dm_annsum)
######En total 9078 segmentos anotados como CDS
## 3265 CDS con el nuevo bed generado con el filtrado

# Count the occurrences of classifications in the Status
# column across the annotation types.
dm_catsum = summarize_categorical(
  annotated_regions = dm_annotated,
  by = c('annot.type', title_str),
  quiet = TRUE)
print(dm_catsum)

dm_catsum2 = summarize_categorical(
  annotated_regions = dm_annotated,
  by = c('annot.type', score),
  quiet = TRUE)
print(dm_catsum2)



#################################################
##############      Normalizar       ############
#################################################
# Normalizamos por el número de segmentos que pertenecen a cada probabilidad
df_calidad <- read.table(file = "intersect_E1_f0.7.bed")
colnames(df_calidad) <- c("chrom", "start", "end", "state", "score")
data <- table(cut(df_calidad$score, breaks = seq(0.65, 1, by = 0.05),
                  labels =paste(c("<0.7", "0.7-0.75", "0.75-0.8", "0.8-0.85",
                                  "0.85-0.9", "0.9-0.95", "0.95-1"))))
# Resultados
# <0.7 0.7-0.75 0.75-0.8 0.8-0.85 0.85-0.9 0.9-0.95   0.95-1 
# 0      862      337      936      841     1275    39757




#################################################
##############         CDS           ############
#################################################
# Nos quedamos con las anotaciones solo de CDS
dm_cds <- dm_catsum2[dm_catsum2$annot.type == "hg19_genes_cds", ]

# Ahora las vamos a agrupar por probabilidad posterior
resum_cds <-cut(dm_cds$Post.P, breaks = seq(0.7, 1, by = 0.05))
table(resum_cds)

# Generamos un barplot de la distribución de la anotación en función de probabilidad posterior
library(RColorBrewer)
pdf("Distribucion_cds.pdf")
coul <- brewer.pal(6, "Set2")
par(mar = c(6, 6, 6, 6), las = 2)
barplot(table(resum_cds), col = coul,
        main = "Distribution of CDS annotation among E1 posterior probability",
        cex.main= 1.3,
        cex.lab = 1, 
        ylab = "Number of segments",
        xlab = "Posterior Probability")
dev.off()

# Normalizamos
# Ahora normalizando
pdf("Porcentaje_cds")
norm_cds <- table(resum_cds)*100 / data[-1]
barplot(norm_cds, col = coul,
        main = "Percentage of segments annotated as promoters",
        cex.main= 1.3,
        cex.lab = 1.3, 
        ylab = "Percentage (%)",
        xlab = "Posterior Probability")
dev.off()

#################################################
##############      Promotores       ############
#################################################
dm_promoter <- dm_catsum2[dm_catsum2$annot.type == "hg19_genes_promoters", ]

resum_promoter <-cut(dm_promoter$Post.P, breaks = seq(0.7, 1, by = 0.05))
table(resum_promoter)
# (0.7,0.75] (0.75,0.8] (0.8,0.85] (0.85,0.9] (0.9,0.95]   (0.95,1] 
# 163        113        247        245        278        456

library(RColorBrewer)
coul <- brewer.pal(6, "Set2")
pdf("Distribucion_promotores.pdf")
par(mar = c(6, 6, 6, 6))
barplot(table(resum_promoter), col = coul,
        main = "Distribution of promoter annotation among E1 posterior probability",
        cex.main= 1.2,
        cex.lab = 1.3, 
        ylab = "Frequency",
        xlab = "Posterior Probability")
dev.off()
# Ahora normalizando
pdf("Porcentaje_promotores")
norm_promoters <- table(resum_promoter)*100 / data[-1]
barplot(norm_promoters, col = coul,
        main = "Porcentaje de segmentos de cada probabilidad posterior que son promotores",
        cex.main= 1.3,
        cex.lab = 1.3, 
        ylab = "Porcentaje",
        xlab = "Posterior Probability")
dev.off()

#View a heatmap of regions occurring in pairs of annotations
annots_order = c(
  'hg19_genes_promoters',
  'hg19_genes_5UTRs',
  'hg19_genes_exons',
  'hg19_genes_introns',
  'hg19_genes_3UTRs',
  'hg19_genes_cds',
  'hg19_enhancers_fantom')
dm_vs_coannotations = plot_coannotations(
  annotated_regions = dm_annotated,
  annotation_order = annots_order,
  axes_label = 'Annotations',
  plot_title = 'Regions in Pairs of Annotations')
print(dm_vs_coannotations)

dm_annotations_plot = plot_annotation(annotated_regions = dm_annotated,
                annotation_order = annots_order,
                plot_title = 'Annotations for E1',
                x_label = 'Annotation type',
                y_label = 'Count')
print(dm_annotations_plot)

# View the number of regions per annotation and include the annotation
# of randomized regions
library("ggplot2")
install.packages("ggthemes") # Install 
library(ggthemes)
pdf("Anotaciones_vs_random.pdf")
dm_annotations_plot_wrandom <-  plot_annotation(
  annotated_regions = dm_annotated,
  annotated_random = dm_random_annotated,
  annotation_order = annots_order,
  plot_title = 'Annotations for E1 (with rndm.)',
  x_label = 'Annotation type',
  y_label = 'Count')

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

dm_annotations_plot_wrandom + theme_hc() +
  scale_fill_hc()

# Otras paletas por probar
  # scale_fill_economist()
  # scale_fill_brewer(palette = "Set2")
ggsave("Anotacion3.pdf")


# dm_annotations_plot_wrandom + theme_stata()+
#   geom_bar(aes(fill = data_type))
  
print(dm_annotations_plot_wrandom)
dev.off()


#CpGIslands
annots_order = c(
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_shelves',
  'hg19_cpg_inter')
dm_annotations_plot_wrandom = plot_annotation(
  annotated_regions = dm_annotated,
  annotated_random = dm_random_annotated,
  annotation_order = annots_order,
  plot_title = 'Annotations for E1 (with rndm.)',
  x_label = 'Annotation type',
  y_label = 'Count')

dm_annotations_plot_wrandom + theme_hc() +
  scale_fill_hc()
ggsave("Anotacion_isla_CpG.pdf")
print(dm_annotations_plot_wrandom)



########################################################################################
########################################################################################
###########TODO ESTO PARA VARIOS ESTADOS O VARIAS CLASES CON ALGUNA CATEGORÍA###########
########################################################################################
########################################################################################
# View the proportions of data classes in knownGene annotations
# The orders for the x-axis labels.
x_order = c(
  'hg19_genes_promoters',
  'hg19_genes_5UTRs',
  'hg19_genes_exons',
  'hg19_genes_introns',
  'hg19_genes_3UTRs',
  'hg19_genes_cds')
dm_vs_kg_cat = plot_categorical(
  annotated_regions = dm_annotated, x='annot.type', fill=title_str,
  x_order = x_order, position='fill',
  legend_title = title_str,
  x_label = 'knownGene Annotations',
  y_label = 'Proportion')
print(dm_vs_kg_cat)

# View the proportions of data classes in CpG annotations
# The orders for the x-axis labels.
x_order = c(
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_shelves',
  'hg19_cpg_inter')
# Make a barplot of the data class where each bar
# is composed of the counts of CpG annotations.
dm_vs_cpg_cat1 = plot_categorical(
  annotated_regions = dm_annotated, x='annot.type', fill=title_str,
  x_order = x_order, position='fill',
  legend_title = title_str,
  x_label = 'CpG Annotations',
  y_label = 'Proportion')
print(dm_vs_cpg_cat1)

# View the counts of CpG annotations in data classes
# The orders for the x-axis labels. This is also a subset
# of the labels (E1, E2, E3...).
x_order = c()
# The orders for the fill labels. Can also use this
# parameter to subset annotation types to fill.
fill_order = c(
  'hg19_genes_promoters',
  'hg19_genes_5UTRs',
  'hg19_genes_exons',
  'hg19_genes_introns',
  'hg19_genes_3UTRs',
  'hg19_genes_cds')
# Make a barplot of the data class where each bar
# is composed of the counts of CpG annotations.
dm_vs_cpg_cat1 = plot_categorical(
  annotated_regions = dm_annotated, x=title_str, fill='annot.type',
  x_order = x_order, fill_order = fill_order, position='stack',
  plot_title = 'Status by Annotation Counts',
  legend_title = 'Annotations',
  x_label = title_str,
  y_label = 'Count')
print(dm_vs_cpg_cat1)
# Make a barplot of the data class where each bar
# is composed of the *proportion* of CpG annotations.
dm_vs_cpg_cat2 = plot_categorical(
  annotated_regions = dm_annotated, x=title_str, fill='annot.type',
  x_order = x_order, fill_order = fill_order, position='fill',
  plot_title = 'Status by CpG Annotation Counts',
  legend_title = 'Annotations',
  x_label = title_str,
  y_label = 'Proportion')
print(dm_vs_cpg_cat2)
