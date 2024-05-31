#!/bin/bash

path=/home/liuxue/data/bioreactor_sludge/segment
cd $path
source activate segment

files=(

0.8

)


for file in "${files[@]}";
    do

for file_1 in `cat list`;
do
$path/segment_find-333.sh --indata /home/liuxue/data/PRJNA605178/beifen  --library_ID $file_1 --taxidDB /home/liuxue/software/prot.accession2taxid  --cor $file  --datatype 2 --nt_noViruses /home/liuxue/database/NT-novirus/nt_noViruses-2023-5-8/nt_noViruses --nt /home/liuxue/software/nt/nt_20221015/nt --nr_virus /home/liuxue/database/NR/nr_virus --nr /home/liuxue/database/NR/nr
done

 cd $path/$file
 cp $path/hebing.R $path/$file
 Rscript hebing.R
 cp $path/contigs_number.R  $path/$file
 Rscript contigs_number.R
 mv total-600.final.wide_data.xlsx total-600-"$file".final.wide_data.xlsx
 mv total-600.final.confidence_table.xlsx total-600-"$file".final.confidence_table.xlsx
 cd $path
 mkdir rdrp_multi10-nordrp10-200
 mv $path/$file rdrp_multi10-nordrp10-200

done;	    


