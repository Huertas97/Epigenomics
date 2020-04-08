#!/bin/bash

path=$1


# Necesita carpeta donde esten los archivos. 
if [ -z "$path" ]
then
      echo "Necesitas introducir el path hasta los archivos"
      echo "Por ejemplo: ./bed_E1_filtrado_0.9/by_chr/"
      exit 123

fi 

# ls $1
files_mon_1=$(ls $1 | grep -E "Monocyte_1_")
files_mon_2=$(ls $1 | grep -E "Monocyte_2_")



samples=$(echo {1..22} X M)

[ -e $path"01_jaccard_score.txt" ] && rm $path"01_jaccard_score.txt"

for i in {1..24}
do 

    file1=$(echo $files_mon_1 | cut -d " " -f$i) 
    file2=$(echo $files_mon_2 | cut -d " " -f$i)

    x=$(echo $file1 | tr "_" "." | cut -d "." -f6)
    echo $x 
    echo $file1
    echo $file2
    jaccard=$(bedtools jaccard -a $path$file1 -b $path$file2)
    echo -e $x"\n$jaccard\n" >> $path"01_jaccard_score.txt"


done 




# bedtools jaccard -a $file1 -b $file2