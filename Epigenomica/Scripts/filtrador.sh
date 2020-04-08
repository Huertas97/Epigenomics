#!/bin/bash

# OBJETIVO: Se encarga de generar un archivo con la intersección de los dos monocitos
# mediante segmentos de 200 bp. Solo se emplearán los segmentos que hayan sido asignados
# al estado deseado (E1, E2 ...)  con una probabiidad posterior deseada.

# EMPLEO: primer argumento el filtro, segundo ele stado deseadp
# ejemplo para coger la interseccion del E1 de ambos monocitos con filtro 0.7
# bash filtrador.sh 0.7 1

filter=$1
n_samples=2
estado=$2

# Filtro por defecto
if [ -z "$filter" ]
then
      echo "Filter not selected. Default value: 0.7"
      filter=0.7
else
    echo "Filter selected: " $filter
fi 

# Numero de muestras por defecto
if [ -z "$n_samples" ]
then
      echo "Number of samples not selected. Default value: 2"
      n_samples=2
else
    echo "Number of samples selected: " $n_samples
fi 

if [ -z "$estado" ]
then
      echo "Number of chromatine state not selected. Default state: 1"
      estado=1
else
    echo "Chromatine state selected: " $estado
fi 

# Generamos carpeta de salida de resultados
mkdir -p "bed_E"$estado"_filtrado_"$filter
mkdir -p "bed_E"$estado"_filtrado_"$filter"/by_chr"

# En cada iteracion cogemos un monocito (monocito 1 y monocito 2)
for i in $(seq 1 $n_samples)
do
    lista_files=$(ls ./RESULTS/Modelo_11_estados/POSTERIOR/ | grep "Monocyte$i")
    # Creamos el archivo summary de num de segmentos por cromosoma comprobando que no existe previamente
    [ -e "bed_E"$estado"_filtrado_"$filter"/by_chr/Monocyte_"$i"chr_summary.txt" ] && rm "bed_E"$estado"_filtrado_"$filter"/by_chr/Monocyte_"$i"chr_summary.txt"

    # Creamos el archivo con la concatenacion de todos los cromosomas por monocito
    [ -e "bed_E"$estado"_filtrado_"$filter"/Monocyte_"$i"_E"$estado"_f"$filter".bed" ] && rm "bed_E"$estado"_filtrado_"$filter"/Monocyte_"$i"_E"$estado"_f"$filter".bed"
    for file in $lista_files
    do 
        
        chr=$(echo $file | grep -oE "chr[0-9MX][0-9]?") # cogemos solo el cromosoma del nombre, siendo num o letra permitiendo un segundo digito con ? (0 o ninguno)
        # echo $chr
        # Generamos el archivo filtrado por cada cromosoma. Genero un archivo para cada cromosoma
        # De las prob posteriores nos quedamos co la tercera fila para adelante quitando el encabezado.
        # Añadimos dos columnas al final que van en ventanas de 200 en 200
        # Nos quedamos luego con las líneas que tienen un valor de probabilidad superior o igual al filtro deseado. El estado es la columna que miramos
        # Luego reconstruyo el documento quedandome solo con el chr, las posiciones y el estado
        tail -n +3 ./RESULTS/Modelo_11_estados/POSTERIOR/$file | awk 'BEGIN {start=0; end=200; line=1} {OFS="\t"; print $0, start,end*line; start=end*line; line+=1}' | awk -v filter=$filter -v estado=$estado '{OFS="\t"; if ($estado >= filter) print $0 }' |  awk -v chr=$chr -v estado=$estado '{OFS="\t"; print chr,$12,$13, "E"estado}' > "bed_E"$estado"_filtrado_"$filter"/by_chr/Monocyte_"$i"_E"$estado"_f"$filter"_"$chr".bed"
        
        # Hacemos un summary por cada cromosoma y lo añadimos a un fichero summary para cada monocito
        echo $file
        wc -l "bed_E"$estado"_filtrado_"$filter"/by_chr/Monocyte_"$i"_E"$estado"_f"$filter"_"$chr".bed" | tr "_" "." | tr " " "." | cut -d "." -f1,12 | tr "." "\t"
        wc -l "bed_E"$estado"_filtrado_"$filter"/by_chr/Monocyte_"$i"_E"$estado"_f"$filter"_"$chr".bed" | tr "_" "." | tr " " "." | cut -d "." -f1,12 | tr "." "\t" >> "bed_E"$estado"_filtrado_"$filter"/by_chr/Monocyte_"$i"chr_summary.txt"
        
        # Generamos el archivo filtrado por monocito uniendo todos los cromosomas. Añado el cromosoma de la iteracion a uno global
        cat "bed_E"$estado"_filtrado_"$filter"/by_chr/Monocyte_"$i"_E"$estado"_f"$filter"_"$chr".bed" >> "bed_E"$estado"_filtrado_"$filter"/Monocyte_"$i"_E"$estado"_f"$filter".bed"
    done

    # Hacemos la media del summary de ambos monocitos
    paste "bed_E"$estado"_filtrado_"$filter"/by_chr/"Monocyte_1chr_summary.txt "bed_E"$estado"_filtrado_"$filter"/by_chr/"Monocyte_2chr_summary.txt | awk 'BEGIN{n=0; OFS="\t"} {n=($1+$3)/2; print n,$2; n=0}' > "bed_E"$estado"_filtrado_"$filter"/by_chr/02_summary_segments.txt"

    # Ahora hacemos la interseccion de los segmentos de todo el genoma de ambos monocitos 
    bedtools intersect  -a "bed_E"$estado"_filtrado_"$filter"/Monocyte_1_E"$estado"_f"$filter".bed" -b "bed_E"$estado"_filtrado_"$filter"/Monocyte_2_E"$estado"_f"$filter".bed" -sorted > "bed_E"$estado"_filtrado_"$filter"/intersect_E"$estado"_f"$filter".bed"
done
