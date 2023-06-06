#!/bin/bash -l
#SBATCH --partition=compute
#SBATCH --mem=50G
#SBATCH --cpus-per-task=50
#SBATCH --mail-user=ayse-oshima@oist.jp
#SBATCH --mail-type=FAIL,END
#SBATCH --time=15

#### Synopsis:
#     1.   Cutadapt trimming off primers, adapter seqs, filtering out N-containing reads, quality cut off etc.
#     2-3. Merging, dereplication, denoising to get ASV table with DADA2 (QIIME plug-in)
# --> 4.   Create database with the makeblastdb the MitoFish seqs 
#     5.   BLAST the ASV table against the makblastdb'd DB made from MitoFish seqs
#     4.   LCA script


helpFunction()
{
   echo ""
   echo "Usage: $0 -o output_directory"
   echo -e "\t-o output file directory for makeblastdb and blastn"
   exit 1 # Exit script after printing help
}

while getopts "o:" opt
do
   case "$opt" in
      o ) blast_dir="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
##for array jobs, use $SLURM_ARRAY_TASK_ID, and submit as "sbatch --array=0-31, --array=1,4,10..."
tid=$SLURM_ARRAY_TASK_ID

# 3. Preparing the BLAST database from the MitoFish fasta file with makeblastdb #######################################################################################################################

# mkdir and cd to the directory for blast work
mkdir -p ${blast_dir%}

cp /flash/RavasiU/Ayse/files_for_makeblastdb/* ${blast_dir%}
cd ${blast_dir%}

# load modules
ml bioinfo-ugrp-modules
ml ncbi-blast/2.10.0+

# modifying files for makeblastdb #############################################
# 1. complete_partial_mitogenomes.fa (downloaded from the MitoFish website)   #
#    http://mitofish.aori.u-tokyo.ac.jp/download/                             #
# 2. nucl_gb.accession2taxid (downloaded from NCBI)                           #
#    https://www.ncbi.nlm.nih.gov/guide/taxonomy/                             #
###############################################################################

# take just the headers of mitogenomes --> mitogenomes_headers
grep "^>" complete_partial_mitogenomes_2023.fa | awk 'sub(/^>/, "")' > mitogenomes_headers

# take just the gb|accession --> mitogenomes_gb_accession
awk -F'|' 'OFS="|" {print $2}' mitogenomes_headers > mitogenomes_gb_accession

# version of complete_partial_mitogenomes.fa where the header is just the accession code --> cleaned_complete_partial_mitogenomes.fa
#        e.g.
##       >gb|KY176037|Prochilodus nigricans (["Characiformes;"] [])
##       ACGTACGATCAGCTAC…
#       becomes:
##      >KY176037
##      ACGTACGATCAGCTAC…
sed 's/>gb|/>/' complete_partial_mitogenomes_2023.fa | sed -e 's/|.*//' > cleaned_complete_partial_mitogenomes.fa

# mapping/meta file that links accession codes to taxid from nucl_gb.accession2taxid
#       Sort both mitogenomes_gb_accession (accession) and nucl_gb.accession2taxid (accession, accession.version, taxid, gi) by gb accession code,
#       then join them by this gb accession code,
#       then output only the gb accession code and the taxID and write into: genbank_taxid
join --nocheck-order -j1 -t $'\t' -o 1.1,2.3   <(sort -t$'\t' -k1 mitogenomes_gb_accession) <(sort -t$'\t' -k1 nucl_gb.accession2taxid) > genbank_taxid


# makeblastdb with cleaned_complete_partial_mitogenomes.fa and genbank_taxid meta file
makeblastdb -in cleaned_complete_partial_mitogenomes.fa -parse_seqids -blastdb_version 5 -taxid_map genbank_taxid -dbtype nucl -title complete_partial_mitogenomes.fa -out ${blast_dir}/complete_partial_mitogenomes/

