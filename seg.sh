#!/bin/bash

path=/home/liuxue/data/PRJDB11882/segment
cd $path

source activate segment


$path/SegFinder.sh --indata /home/liuxue/data/PRJDB11882/data --taxidDB /home/liuxue/software/prot.accession2taxid --nt_noViruses /home/liuxue/database/NT-novirus/nt_noViruses-2023-5-8/nt_noViruses --nt /home/liuxue/software/nt/nt_20221015/nt  --thread 20 --datatype 2  --method salmon --preprocess true  --assemble megahit  --nr /home/liuxue/database/NR/nr --only_rdrp_find 1


#for file in `cat list`;
#do
#$path/SegFinder.sh --indata /home/liuxue/data/bioreactor_sludge/data  --taxidDB /home/liuxue/software/prot.accession2taxid --nt_noViruses /home/liuxue/database/NT-novirus/nt_noViruses-2023-5-8/nt_noViruses --nt /home/liuxue/software/nt/nt_20221015/nt  --thread 20  --rm_length 600 --datatype 2 --cor 0.79 --library_ID $file --method salmon  --nr /home/liuxue/database/NR/nr 

#done
conda deactivate
