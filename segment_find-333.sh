usage(){
cat <<'EOF'
Usage
 [-o]... the directory to output the results; default ./
 [--indata]... the location of the raw data
 [--incontig]... the contig you want to search
 [--thread]... default 10
 [--cor]... correlation coefficient;default 0.8
 [--nt_noViruses]... the location of nt_noViruses database
 [--nt]... the location of nt database
 [--nr]... the location of nr database
 [--nr_virus]... the location of nr_virus database
 [--method]... the method to quantify the transcript abundances,salmon or RSEM,default salmon
 [--datatype]... the type of input data，single(input 1) or double(input 2)
 [--taxidDB]... the location of prot.accession2taxid database
 [--rm_length]... the contigs whose length less than this value will be removed,default 600
 [--min_multi]... if the length of the contigs and their re_assembled cotigs are all less than this value, the clusters they are in will be removed,default 20
 [--library_ID]... the library you want to search,default the library in which the abundance of the contig you input is max
 [--preprocess]... whether to preprocess the raw fq data,you can input 'true' or 'false',default false
 [--assemble]... the tool you choose to assemble the raw reads,megahit or spades,default spades
 [--only_rdrp_find]...1 or 0, 1 means only find rna virus RdRPs without any other analysis, default 0
 [--min_TPM]...if there exist the contig whose TPM is less than this value, the cluster it is in will be removed,default 200
 [--help]...
 [--version]...
EOF
}

parameters=`getopt  -o o: --long indata:,incontig:,thread:,cor:,datatype:,nt:,nr:,nr_virus:,method:,only_rdrp_find:,min_multi:,min_TPM:,assemble:,preprocess:,library_ID:,rm_length:,nt_noViruses:,taxidDB:,help,version -n "$0" -- "$@"`
[ $? -ne 0 ] && { echo "Try '$0 --help' for more information."; exit 1; }

eval set -- "$parameters"
out_loc_flag=0
indata_flag=0
incontig_flag=0
datatype_flag=0
nt_flag=0
nr_flag=0
nr_virus_flag=0
nt_noViruses_flag=0
taxidDB_flag=0
thread=10
cor=0.8
rm_length=600
min_multi=10
min_TPM=100
library_ID_flag=0
quantify_method=salmon
preprocess=false
assemble_method=spades
only_rdrp_find=0
while true;do
    case "$1" in
        --indata)  indata_flag=1; rawData_loc=$2;  shift 2;;
        --incontig) incontig_flag=1; contig=$2; shift 2 ;;
        --thread) thread=$2; shift 2 ;;
        --cor) cor=$2; shift 2 ;;
                --rm_length) rm_length=$2; shift 2 ;;
                --min_multi) min_multi=$2; shift 2 ;;
                --min_TPM) min_TPM=$2; shift 2 ;;
        --nt)  nt_flag=1; nt_loc=$2;  shift 2;;
                --nr)  nr_flag=1; nr_loc=$2;  shift 2;;
                --nr_virus)  nr_virus_flag=1; nr_virus_loc=$2;  shift 2;;
        --nt_noViruses)  nt_noViruses_flag=1; nt_noViruses_loc=$2;  shift 2;;
                --library_ID)  library_ID_flag=1; library_ID=$2;  shift 2;;
                --datatype)  datatype_flag=1; datatype=$2;  shift 2;;
                --method)   quantify_method=$2;  shift 2;;
                --assemble)   assemble_method=$2;  shift 2;;
                --only_rdrp_find)   only_rdrp_find=$2;  shift 2;;
                --preprocess)   preprocess=$2;  shift 2;;
        --taxidDB)  taxidDB_flag=1; taxidDB_loc=$2;  shift 2;;
        -o) out_loc_flag=1; out_loc=$2; shift 2 ;;
        --version) echo "$0 version V1.0"; exit ;;
        --help) usage;exit ;;
        --)
            shift
                        if [[ $preprocess != true && $preprocess != false ]]; then echo "please input true or false!!! --preprocess";  exit 1; fi;
            if [ $indata_flag -eq 0 ]; then echo "please input the location of the raw data!!! --indata"; exit 1; fi;
            if [[ $incontig_flag -eq 0 && $library_ID_flag -eq 0 && $only_rdrp_find -eq 0 ]]; then echo "please input a contig name or a library_ID!!! --incontig or --library_ID";  exit 1; fi;
                        if [[ $nt_flag -eq 0 && $only_rdrp_find -eq 0 ]]; then echo "please input the location of nt!!! --nt";  exit 1; fi;
                        if [[ $nr_flag -eq 0 && $preprocess == true ]]; then echo "please input the location of nr!!! --nr";  exit 1; fi;
                        if [[ $nr_virus_flag -eq 0 && $preprocess == true ]]; then echo "please input the location of nr_virus!!! --nr_virus";  exit 1; fi;
                        if [[ $nt_noViruses_flag -eq 0 && $only_rdrp_find -eq 0 ]]; then echo "please input the location of nt_noViruses database!!! --nt_noViruses";  exit 1; fi;
                        if [ $datatype_flag -eq 0 ]; then echo "please input the type of input data!!! --datatype";  exit 1; fi;
                        if [ $taxidDB_flag -eq 0 ]; then echo "please input the location of prot.accession2taxid database!!! --taxidDB";  exit 1; fi;
            if [ $out_loc_flag -eq 0 ]; then echo "the output files are all in the current directory"; fi
            break ;;
        *) usage;exit 1;;
    esac
done
if [[ $datatype -ne 1 && $datatype -ne 2 ]]; then echo 'please re_input the type of input data, 1 or 2,for 1 means single type,2 means double';  exit 1; fi
if [ $out_loc_flag -eq 0 ]; then megahit=megahit;  nr=nr; rdrp=rdrp; network=network; fi
if [ $out_loc_flag -eq 1 ]; then megahit=$out_loc/megahit; nr=$out_loc/nr; rdrp=$out_loc/rdrp; network=$out_loc/network; fi

chmod +x align_and_estimate_abundance.pl
chmod +x ORFfinder
#################################
present_loc=`pwd`
megahit_loc=$present_loc/$megahit
path=$rawData_loc




	rm -rf ${present_loc}/${library_ID}/megahit/${library_ID}.final.confidence_table.xlsx
rm -rf ${present_loc}/${library_ID}/megahit/${library_ID}.pre.confidence_table.xlsx
rm -rf ${present_loc}/${library_ID}/megahit/${library_ID}.network_group_fr.pdf
rm -rf ${present_loc}/${library_ID}/network/*
cp Cor_contigs_extract.R ${library_ID}/$megahit/
cd ${library_ID}/$megahit

awk '/^>/ {printf("%s\t", substr($0,2)); getline; print length($0)}' "${library_ID}".re.fasta  | awk '{match($0, /len.*[0-9]*/); str=substr($0, RSTART, RLENGTH); match(str, /[0-9]+/); a=substr(str, RSTART, RLENGTH);if (a == $2) print $0"\tfalse"; else print $0"\ttrue"}' > re.fasta_length.txt
cat re.fasta_length.txt |uniq > re.fasta_length.txt-1
mv  re.fasta_length.txt-1 re.fasta_length.txt
Rscript Cor_contigs_extract.R  ${cor} ${library_ID} 200 10 10


cd ${present_loc}
#seqtk subseq ${present_loc}/megahit/${library_ID}.megahit.fa  ${present_loc}/megahit/Cor_contigs.txt  > ${present_loc}/megahit/${library_ID}.re_contigs.fasta
mkdir ${cor}
#awk '{if($0~/^>/)a=$0; else{system("echo \\"a"\"\n\""$0a)}}' ${present_loc}/megahit/${library_ID}.re_contigs.fasta
#mv ${present_loc}/${library_ID}* $megahit
cp ${present_loc}/${library_ID}/megahit/${library_ID}.network_group_fr.pdf ${present_loc}/${library_ID}/network/
cp ${present_loc}/${library_ID}/megahit/${library_ID}.final.confidence_table.xlsx ${present_loc}/${library_ID}/network/
cp ${present_loc}/${library_ID}/megahit/${library_ID}.pre.confidence_table.xlsx ${present_loc}/${library_ID}/network/
cp ${present_loc}/${library_ID}/network/${library_ID}* ${present_loc}/${cor}
