#!/bin/bash

path=/home/liuxue/data/plants+/segment
cd $path

source activate segment
for file in `cat list`;
do
/home/liuxue/data/PRJNA605178/segment/segment_find-333.sh --indata /home/liuxue/data/PRJNA605178/data  --library_ID $file --taxidDB /home/liuxue/software/prot.accession2taxid  --cor 0.8  --datatype 2 --nt_noViruses /home/liuxue/database/NT-novirus/nt_noViruses-2023-5-8/nt_noViruses --nt /home/liuxue/software/nt/nt_20221015/nt --nr_virus /home/liuxue/database/NR/nr_virus --nr /home/liuxue/database/NR/nr

done
conda deactivate
