#!/bin/bash

#Necesario el fichero con el genoma.
#Para descargarlo:
#     wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz 
#     gzip -d hg19.fa.gz

file_in=$1
file_out=$2

bedtools getfasta -fi hg19.fa -bed $file_in -fo $file_out 

#Ejemplo
# ./getFastaFromBed_MEME.sh bed_E1_filtrado_0.7/intersect_E1_f0.7.bed bed_E1_filtrado_0.7/intersect_E1_f0.7.fasta