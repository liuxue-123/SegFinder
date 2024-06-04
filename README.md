# SegFinder
## segmented virus finder workflow
![](https://github.com/liuxue-123/SegFinder/blob/main/flow/workflow.png)

SegFinder detection pipeline. a, Schematic overview of the discovery of RdRP for RNA viruses. The inputs are fastq files for multiple meta-transcriptome libraries. rRNA, ribosomal RNA; NR, Non-Redundant Protein Sequence Database; NT, Nucleotide Sequence Database. b, The processing pipeline of correlation calculation. L, library; C, contig; c, Schematic illustration of filtering of segmented RNA virus clusters. Cor, correlation; TPM, Transcripts Per Kilobase of exon model per Million mapped reads;

## 1.Installation
  ### 1.1 Install conda and SegFinder dependencies

```
conda env create -f environment.yml
conda activate SegFinder
```

### 1.2 Downloading and configuring the database

  #### 1.2.1 prot.accession2taxid

```
mkdir Seg_DB
cd Seg_DB
mkdir accession2taxid
cd accession2taxid
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5

#Check for the file integrity
md5sum -c prot.accession2taxid.gz.md5

#Unzip the files and onfiguration
gunzip -c prot.accession2taxid.gz > Seg_DB/accession2taxid/prot.accession2taxid
```


  #### [1.2.2 NCBI Non-Redundant Protein Database (NR)](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)
```
cd Seg_DB
mkdir nr
cd nr
wget -t 0 -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
```

  #### [1.2.3 NCBI Nucleotide Sequence Database (NT)](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)
```
cd Seg_DB
mkdir nt
cd nt
wget -t 0 -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
```
  #### 1.2.4 Virus-free non-redundant nucleotide (virus-free nt)
```
#download nucl_gb.accession2taxid
cd Seg_DB/accession2taxid
wget -t 0 -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip -c nucl_gb.accession2taxid.gz > Seg_DB/accession2taxid/nucl_gb.accession2taxid

#download taxdump
mdkir Seg_DB/taxdump
wget -t 0 -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
cd Seg_DB/taxdump
tar -zxvf taxdump.tar.gz

#download ref_viruses_rep_genomes
mkdir Seg_DB/ref_viruses_rep_genomes
cd Seg_DB/ref_viruses_rep_genomes
wget -t 0 -c https://figshare.com/ndownloader/files/46795402
tar -zxvf ref_viruses_rep_genomes.tar.gz

#Handling nt database
python3 process_sequences.py --input Seg_DB/nt/ --out Seg_DB/nt/ --threads 40  --nucl_gb_accession2taxid_path Seg_DB/accession2taxid/nucl_gb.accession2taxid --taxdump_path Seg_DB/taxdump --ref_viruses_path Seg_DB/ref_viruses_rep_genomes
```
Note: --input:nt database location 

## 2.Usage

#### Step 1: discovery of RdRP for RNA viruses 

./SegFinder.sh [option] --help  

```
./SegFinder.sh --indata PATH/data \
               --taxidDB Seg_DB_PATH/prot.accession2taxid \
               --nt_noViruses PATH/nt_noViruses \
               --nt Seg_DB_PATH/nt \
               --thread 20 \
               --datatype 2 \
               --method salmon \
               --preprocess true  \
               --assemble megahit  \
               --nr Seg_DB_PATH/nr \
               --only_rdrp_find 1
```
Note:file_list.txt contains the prefix name of the file;The file name must be prefixed_1/2.fq.gz format
#### Step 2: segmented RNA virus finder 
```
./SegFinder.sh --indata PATH/data \
               --taxidDB Seg_DB_PATH/prot.accession2taxid \
               --nt_noViruses Seg_DB_PATH/nt_noViruses \
               --nt Seg_DB_PATH/nt  \
               --thread 20 \
               --rm_length 600 \
               --datatype 2 \
               --cor 0.8 \
               --library_ID $file \
               --method salmon  \
               --nr Seg_DB_PATH/nr
Note:list.txt contains the prefix name of the file;The file name must be prefixed_1/2.fq.gz format
```  
