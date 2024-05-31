# SegFinder
## segmented virus finder workflow
![](https://github.com/liuxue-123/SegFinder/blob/main/flow/workflow.png)

SegFinder detection pipeline. a, Schematic overview of the discovery of RdRP for RNA viruses. The inputs are fastq files for multiple meta-transcriptome libraries. rRNA, ribosomal RNA; NR, Non-Redundant Protein Sequence Database; NT, Nucleotide Sequence Database. b, The processing pipeline of correlation calculation. L, library; C, contig; c, Schematic illustration of filtering of segmented RNA virus clusters. Cor, correlation; TPM, Transcripts Per Kilobase of exon model per Million mapped reads;

## 1.Installation
Step 1: Install conda and SegFinder dependencies

```conda env create -f environment.yml```

Step 2: Downloading and configuring the database


### PROT_ACC2TAXID
#### Download the `PROT_ACC2TAXID` file
``` wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz ```
```wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5```

#### Check for the file integrity
```md5sum -c prot.accession2taxid.gz.md5```

#### Unzip the files and onfiguration
```gunzip -c prot.accession2taxid.gz > PHYRVM_DB_PATH/accession2taxid/prot.accession2taxid```

#### NCBI Non-Redundant Protein Database (NR)
#### NCBI Nucleotide Sequence Database (NT)
#### Virus-free non-redundant nucleotide (virus-free nt)


## Usage
./SegFinder.sh --indata /home/liuxue/data/PRJDB11882/data --taxidDB /home/liuxue/software/prot.accession2taxid --nt_noViruses /home/liuxue/database/NT-novirus/nt_noViruses-2023-5-8/nt_noViruses --nt /home/liuxue/software/nt/nt_20221015/nt  --thread 20 --datatype 2  --method salmon --preprocess true  --assemble megahit  --nr /home/liuxue/database/NR/nr --only_rdrp_find 1

./SegFinder.sh --indata /home/liuxue/data/bioreactor_sludge/data  --taxidDB /home/liuxue/software/prot.accession2taxid --nt_noViruses /home/liuxue/database/NT-novirus/nt_noViruses-2023-5-8/nt_noViruses --nt /home/liuxue/software/nt/nt_20221015/nt  --thread 20  --rm_length 600 --datatype 2 --cor 0.79 --library_ID $file --method salmon  --nr /home/liuxue/database/NR/nr 
