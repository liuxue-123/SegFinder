import os
import subprocess
import argparse

def run_command(command):
    subprocess.run(command, shell=True, check=True)

def process_sequences(input_path, out_path, num_threads, accession2taxid_path, taxdump_path, ref_viruses_path):
    # 定义路径
    nt_fasta = os.path.join(input_path, "nt")
    nt_fasta_id = os.path.join(out_path, "nt.fasta.id")
    nt_fasta_id_noVirus = os.path.join(out_path, "nt.fasta.id_noVirus")
    nt_fasta_id_noVirus_accNum = os.path.join(out_path, "nt.fasta.id_noVirus_accNum")
    nt_fasta_id_noVirus_accNum_taxid = os.path.join(out_path, "nt.fasta.id_noVirus_accNum_taxid")
    nt_fasta_id_noVirus_accNum_taxid_results = os.path.join(out_path, "nt.fasta.id_noVirus_accNum_taxid_results")
    nt_noVirus_fasta = os.path.join(out_path, "nt_noVirus.fasta")
    id2_txt = os.path.join(out_path, "id2.txt")
    id3_txt = os.path.join(out_path, "id3.txt")
    id3_accessionID = os.path.join(out_path, "id3.accessionID")
    id3_accessionID_fasta = os.path.join(out_path, "id3.accessionID.fasta")
    id3_accessionID_fasta_csv = os.path.join(out_path, "id3.accessionID.fasta.csv")
    id3_virus_accession = os.path.join(out_path, "id3.virus.accession")
    id3_virus_accessionID = os.path.join(out_path, "id3.virus.accessionID")
    id3_diff_txt = os.path.join(out_path, "id3-diff.txt")
    id3_novirus_accessionID_fasta = os.path.join(out_path, "id3-novirus.accessionID.fasta")
    nt_noVirus_update_fasta = os.path.join(out_path, "nt_noVirus-update.fasta")
    nt_noViruses = os.path.join(out_path, "nt_noViruses")

    # 使用 seqkit 提取序列名称
    run_command(f"seqkit seq -j {num_threads} -n {nt_fasta} > {nt_fasta_id}")


    # 使用 grep 过滤掉包含 "virus" 的行
    run_command(f"grep -v -i 'virus' {nt_fasta_id} > {nt_fasta_id_noVirus}")

    # 进一步过滤掉其他病毒相关的关键词
    keywords = ["viruses", "riboviria", "phage", "viridae", "proviral"]
    temp_file = nt_fasta_id_noVirus + "_tmp"
    run_command(f"cp {nt_fasta_id_noVirus} {temp_file}")
    for keyword in keywords:
        run_command(f"grep -v -i '{keyword}' {temp_file} > {temp_file}_filtered")
        run_command(f"mv {temp_file}_filtered {temp_file}")
    run_command(f"mv {temp_file} {nt_fasta_id_noVirus}")

    # 提取序列号并存储在 nt.fasta.id_noVirus_accNum
    run_command(f"sed 's/ /\\t/g' {nt_fasta_id_noVirus} | cut -f1 > {nt_fasta_id_noVirus_accNum}")

    # 使用 grep 进行过滤
    run_command(f"grep -F -f {nt_fasta_id_noVirus_accNum} -w {accession2taxid_path} > {nt_fasta_id_noVirus_accNum_taxid}")

    # 使用 taxonkit 获取分类信息
    run_command(f"cut -f3 {nt_fasta_id_noVirus_accNum_taxid} | sort -u | taxonkit --data-dir {taxdump_path} lineage | awk '$2>0' | taxonkit --data-dir {taxdump_path} reformat -f '{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}' -F | cut -f1,3- > {nt_fasta_id_noVirus_accNum_taxid_results}")

    # 过滤病毒相关的结果
    run_command(f"grep -v -i 'viruses' {nt_fasta_id_noVirus_accNum_taxid_results} | cut -f1 > {os.path.join(out_path, 'nt.fasta.id_noVirus_accNum_taxid_results_noVirusID')}")
    run_command(f"grep -F -f {os.path.join(out_path, 'nt.fasta.id_noVirus_accNum_taxid_results_noVirusID')} -w {nt_fasta_id_noVirus_accNum_taxid} | cut -f2 > {os.path.join(out_path, 'nt.fasta.id_noVirus_accNum_taxid_results_noVirusAccession')}")
 # 提取非病毒序列
    run_command(f"seqtk subseq {nt_fasta} {os.path.join(out_path, 'nt.fasta.id_noVirus_accNum_taxid_results_noVirusAccessionID')} > {nt_noVirus_fasta}")
    run_command(f"cut -f2 {nt_fasta_id_noVirus_accNum_taxid} | sort | uniq > {id2_txt}")

    # 执行 diff.py 的功能
    with open(nt_fasta_id_noVirus_accNum, "r") as file1, open(id2_txt, "r") as file2:
        diff = set(file1).difference(file2)
    with open(id3_txt, "w") as output:
        for line in diff:
            output.write(line)

    # 进一步过滤序列并进行 BLAST 比对
    run_command(f"grep -F -f {id3_txt} -w {nt_fasta_id_noVirus} > {id3_accessionID}")
    run_command(f"seqtk subseq {nt_fasta} {id3_accessionID} > {id3_accessionID_fasta}")
    run_command(f"blastn -query {id3_accessionID_fasta} -db {ref_viruses_path} -out {id3_accessionID_fasta_csv} -evalue 1E-3 -outfmt '6 qseqid sacc staxid salltitles pident evalue' -max_target_seqs 1 -num_threads {num_threads}")
    run_command(f"cut -f1 {id3_accessionID_fasta_csv} | sort | uniq > {id3_virus_accession}")
    run_command(f"grep -F -f {id3_virus_accession} -w {id3_accessionID} > {id3_virus_accessionID}")

    # 执行 diff-1.py 的功能
    with open(id3_accessionID, "r") as file1, open(id3_virus_accessionID, "r") as file2:
        diff = set(file1).difference(file2)
    with open(id3_diff_txt, "w") as output:
        for line in diff:
            output.write(line)

    # 提取最终的非病毒序列
    run_command(f"seqtk subseq {id3_accessionID_fasta} {id3_diff_txt} > {id3_novirus_accessionID_fasta}")

    # 合并结果并创建新的 BLAST 数据库
    run_command(f"cat {nt_noVirus_fasta} {id3_novirus_accessionID_fasta} > {nt_noVirus_update_fasta}")
    run_command(f"makeblastdb -in {nt_noVirus_update_fasta} -dbtype nucl -out {nt_noViruses} -parse_seqids")


    # 删除中间文件
    intermediate_files = [
        nt_fasta_id, nt_fasta_id_noVirus, nt_fasta_id_noVirus_accNum, nt_fasta_id_noVirus_accNum_taxid,
        nt_fasta_id_noVirus_accNum_taxid_results, nt_noVirus_fasta, id2_txt, id3_txt,
        id3_accessionID, id3_accessionID_fasta, id3_accessionID_fasta_csv, id3_virus_accession,
        id3_virus_accessionID, id3_diff_txt, id3_novirus_accessionID_fasta,
        os.path.join(out_path, 'nt.fasta.id_noVirus_accNum_taxid_results_noVirusAccessionID'),
        os.path.join(out_path, 'nt.fasta.id_noVirus_accNum_taxid_results_noVirusID'),
        os.path.join(out_path, 'nt.fasta.id_noVirus_accNum_taxid_results_noVirusAccession')
    ]
    for file in intermediate_files:
        if os.path.exists(file):
            os.remove(file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process sequences and filter out viral sequences.")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads to use.")
    parser.add_argument("--input", type=str, default=".", help="Input path for nt files. Default is current directory.")
    parser.add_argument("--out", type=str, default=".", help="Output path for results. Default is current directory.")
    parser.add_argument("--nucl_gb_accession2taxid_path", type=str, required=True, help="Path to nucl_gb.accession2taxid file.")
    parser.add_argument("--taxdump_path", type=str, required=True, help="Path to taxdump directory.")
    parser.add_argument("--ref_viruses_path", type=str, required=True, help="Path to ref_viruses_rep_genomes database.")

    args = parser.parse_args()
    process_sequences(args.input, args.out, args.threads, args.nucl_gb_accession2taxid_path, args.taxdump_path, args.ref_viruses_path)

