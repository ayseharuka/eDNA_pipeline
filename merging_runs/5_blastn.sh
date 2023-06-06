#!/bin/bash -l
#SBATCH --partition=compute
#SBATCH --mem=50G
#SBATCH --cpus-per-task=50
#SBATCH --mail-user=ayse-oshima@oist.jp
#SBATCH --mail-type=FAIL,END
#SBATCH --time=45

#### Synopsis:
#     1.   Cutadapt trimming off primers, adapter seqs, filtering out N-containing reads, quality cut off etc.
#     2-3. Merging, dereplication, denoising to get ASV table with DADA2 (QIIME plug-in)
#     4.   Create database with the makeblastdb the MitoFish seqs
# --> 5.   BLAST the ASV table against the makblastdb'd DB made from MitoFish seqs
#     6.   LCA script

helpFunction()
{
   echo ""
   echo "Usage: $0 -o output_directory -q query_asv_table.fasta -p perc_id -e e_value"
   echo -e "\t-o output file directory for both  makeblastdb and blastn"
   echo -e "\t-q the asv table that will be given to blastn as a query"
   echo -e "\t-p percent identity cut off"
   echo -e "\t-e expect value for saving hits - lower with more significant matches"
   exit 1 # Exit script after printing help
}

while getopts "o:q:p:e:" opt
do
   case "$opt" in
      o ) blast_dir="$OPTARG" ;;
      q ) query_asv_table="$OPTARG" ;;
      p ) perc_id="$OPTARG" ;;
      e ) e_value="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# 5. BLAST of ASV's from ASV table against the MitoFish db with blastn ######################################################################################################################################

# mkdir and cd to the directory for blast work
path2db=${blast_dir}/complete_partial_mitogenomes
cd ${blast_dir%}

# load modules
ml bioinfo-ugrp-modules
ml ncbi-blast/2.10.0+

blastn -task blastn -db ${path2db} -query ${query_asv_table} -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -out ${blast_dir}/blast_perc_id${perc_id}_eval${e_value}.tab -perc_identity ${perc_id} -evalue ${e_value}

# pid = 95
# e-value has been kept at 1e-2 so far --> wonder if i need to give it as a string??


