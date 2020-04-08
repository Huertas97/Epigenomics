library('annotatr')
title_str <- 'Status'

#Load BED file
#Fichero de prueba solo con 2000 lineas del fichero de segmentos del Monocito 1.
dm_regions = read_regions(con = 'bed_E1_filtrado_0.7/intersect_E1_f0.7.bed', genome = 'hg19', format = 'bed',
                          rename_name = title_str)

# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
annots = c('hg19_genes_cds', 'hg19_basicgenes', 'hg19_cpgs', 'hg19_genes_intergenic',
           'hg19_enhancers_fantom')
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

# Count the occurrences of classifications in the Status
# column across the annotation types.
dm_catsum = summarize_categorical(
  annotated_regions = dm_annotated,
  by = c('annot.type', title_str),
  quiet = TRUE)
print(dm_catsum)

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
dm_annotations_plot_wrandom = plot_annotation(
  annotated_regions = dm_annotated,
  annotated_random = dm_random_annotated,
  annotation_order = annots_order,
  plot_title = 'Annotations for E1 (with rndm.)',
  x_label = 'Annotation type',
  y_label = 'Count')
print(dm_annotations_plot_wrandom)

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
