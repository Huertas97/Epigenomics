#!/bin/bash
mkdir -p "archivos_bed" # no error if existing, make parent directories as needed

estado=$1

# head ./RESULTS/Modelo_11_estados/Monocyte1_11_Master_11_segments.bed
# -wb to report the original “B” feature when an overlap is found
bedtools intersect -wb -a ./RESULTS/Modelo_11_estados/Monocyte1_11_Master_11_segments.bed -b ./RESULTS/Modelo_11_estados/Monocyte2_11_Master_11_segments.bed -sorted  | awk -v var=$estado '{OFS="\t"; if ($4 == $8 && $4 == var) print $1, $2, $3, $4}' > ./archivos_bed/$estado.bed

lines=$(cat ./archivos_bed/$estado.bed | awk '{ $(NF+1)=$3-$2; $(NF+2)=$NF/200; print $0}' | awk 'BEGIN {n=0} {n+=$NF} END {print n}')
echo "Numero de segmentos de 200 bp:" $lines

