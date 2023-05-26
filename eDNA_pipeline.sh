#!/bin/bash -l
#SBATCH --partition=compute
#SBATCH --mem=20G
#SBATCH --cpus-per-task=50
#SBATCH --mail-user=user-name@email.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --time=45

#### Synopsis:
# 1. Cutadapt trimming off primers, adapter seqs, filtering out N-containing reads, quality cut off etc.
# 2. Merging, dereplication, denoising to get ASV table with DADA2 (QIIME plug-in)
# 3. BLAST the ASV table against the makblastdb'd DB made from MitoFish seqs
# 4. LCA script

# Load the necessary modules (loaded in the order of their usage):
ml bioinfo-ugrp-modules DebianMed
ml cutadapt/3.2-2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
##for array jobs, use $SLURM_ARRAY_TASK_ID, and submit as "sbatch --array=0-31, --array=1,4,10..."
tid=$SLURM_ARRAY_TASK_ID

# the delimiter I will use to awk for the basename
FOR=_S

work_dir=path/to/your/working/directory

mkdir -p ${work_dir%}
cd ${work_dir%}

dir=path/to/your/fastq/files

# minimum reads
min_reads=10
# quality threshold
qv=10
## output prefix
out="12S"

# make an output directory for cutadapt
mkdir ${work_dir}/cutadapt_out

####################################################################### Start of the pipeline  #################################################################################################################

# 1. CutAdapt ################################################################################################################################################################################################# 
# create the cutadapt sanity check file with the header (id     r1_before       r1_after        r2_before       r2_after)
echo -e "id\tr1_before\tr1_after\tr2_before\tr2_after" > cutadapt_before_after.txt

#### Loop for grouping R1 and R2 files.

# gets sample base names (id): A1, A2, A3, A-Field etc.
for id in `ls ${dir} |awk -F ${FOR%} '{print $1}' |sort |uniq`
do

# grep for the base names (id) and R1 and R2 to get list of R1 and R2 fastq file names
R1=`ls ${dir} |grep $id |grep R1`
R2=`ls ${dir} |grep $id |grep R2`

## Cutadapt (order matters, will follow the cutadapt documentation)
## (P5 flowcell adapter)-(IDX)-(overhang adapter for 2nd PCR)-(MiFish primer)

## Quality trimming is done before any adapter trimming.
## Adapter trimming itself does not appear in that list and is done after quality trimming and before length trimming (--length/-l).

## explaining cutadapt arguments:
## Processing
#### 1. -u length (+) from the beginning, (-) from the end
#### 2. -q quality cut off
#### 3. -a -g -b adapter/primer
#### 4. -l length (+) from the end, (-) from the beginning
## Filtering of processed reads: too short or untrimmed reads
#### 5. -m discard reads shorter than this -M discard reads longer than this
#### 6. --max-n discard reads with count of N (ambiguous base) greater than the given value. --> set to 0 to remove reads with N>0.

## 1. Filters out sequences with ambiguous (N) bases (DADA2 requires this)
## 2. Removes primers and adapter sequences, as well as filter out qv< values.

# Processing (remove 6 bases form beginning, quality cutoff at 10 from 3 and 5-prime ends, remove the primers, length trim from end for the reverse reads?)

# command for processing and filtering
cutadapt --cores ${OMP_NUM_THREADS} -u 6 -q ${qv},${qv} -g GTCGGTAAAACTCGTGCCAGC -g GCCGGTAAAACTCGTGCCAGC -g RGTTGGTAAATCTCGTGCCAGC -a CAAACTGGGATTAGATACCCCACTATG -a CAAACGGGGATTAGACACCCTCCTATG -a CAAACTAGGATTAGATACCCCACTATGC -G CATAGTGGGGTATCTAATCCCAGTTTG -G CATAGGAGGGTGTCTAATCCCCGTTTG -G GCATAGTGGGGTATCTAATCCTAGTTTG -A GCTGGCACGAGTTTTACCGAC -A GCTGGCACGAGTTTTACCGGC -A GCTGGCACGAGATTTACCAACY -n 2 -m 50 --max-n 0 -o ${work_dir}/cutadapt_out/${id}_trimmed_R1.fastq -p ${work_dir}/cutadapt_out/${id}_trimmed_R2.fastq ${dir}/${R1} ${dir}/${R2}

echo OK1-Cutadapt-Processing

# append the the number of reads before and after the cutadapt to the sanity check file

# input fastq files:
# R1
# R2

# output fq files:
#${work_dir}/trimmed/${id}_trimmed_R1.fq
#${work_dir}/trimmed/${id}_trimmed_R2.fq

before_r1=`zgrep "@" ${dir}/${R1} | wc -l`
before_r2=`zgrep "@" ${dir}/${R2} | wc -l`
after_r1=`grep "@" ${work_dir}/cutadapt_out/${id}_trimmed_R1.fastq | wc -l`
after_r2=`grep "@" ${work_dir}/cutadapt_out/${id}_trimmed_R2.fastq | wc -l`
echo -e $id$'\t'$before_r1$'\t'$after_r1$'\t'$before_r2$'\t'$after_r2 >> cutadapt_before_after.txt

done

# 2. Running DADA2 through the R scripts dada2.r #####################################################################################################################################################

# exporting paths to R and package libraries.
export PATH="/apps/unit/RavasiU/ayse/R-4.2.1/bin:$PATH"
export PATH="/home/a/ayse-oshima/R/x86_64-redhat-linux-gnu-library/4.1:$PATH"

# call/source the R script
/apps/unit/RavasiU/ayse/R-4.2.1/bin/Rscript --vanilla dada2.R

# output files from dada2 will be in: /flash/RavasiU/Ayse/eDNA_workspace/dada2_out
# the files used for downstream steps are:
## asv_map
## asv_table.fasta
## asv_table_mapped.tab
## filtering_results.csv
## seqtab_nochim.csv

# 3. Preparing the BLAST database from the MitoFish fasta file with makeblastdb #######################################################################################################################

# mkdir and cd to the directory for blast work
blast_dir=/path/for/blast/output
mkdir -p ${blast_dir%}

# make sure taxdb files downloaded from 
# https://ftp.ncbi.nlm.nih.gov/blast/db/
# and are in the blast output directory
# I have it in my personal directory from which I copy paste each time.
# wget would work too? but upon updates, the link may not work anymore.

cp /path/for/taxdb/* ${blast_dir%}
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

complete_partial_mitogenomes=/path/for/complete_partial_mitogenomes.fa
nucl_gb_accession2taxid=/path/for/nucl_gb.accession2taxid

# take just the headers of mitogenomes --> mitogenomes_headers
grep "^>" ${complete_partial_mitogenomes} | awk 'sub(/^>/, "")' > mitogenomes_headers

# take just the gb|accession --> mitogenomes_gb_accession
awk -F'|' 'OFS="|" {print $2}' mitogenomes_headers > mitogenomes_gb_accession

# version of complete_partial_mitogenomes.fa where the header is just the accession code --> cleaned_complete_partial_mitogenomes.fa
#        e.g.
##       >gb|KY176037|Prochilodus nigricans (["Characiformes;"] [])
##       ACGTACGATCAGCTAC…
#       becomes:
##      >KY176037
##      ACGTACGATCAGCTAC…
sed 's/>gb|/>/' ${complete_partial_mitogenomes} | sed -e 's/|.*//' > cleaned_complete_partial_mitogenomes.fa

# mapping/meta file that links accession codes to taxid from nucl_gb.accession2taxid
#       Sort both mitogenomes_gb_accession (accession) and nucl_gb.accession2taxid (accession, accession.version, taxid, gi) by gb accession code,
#       then join them by this gb accession code,
#       then output only the gb accession code and the taxID and write into: genbank_taxid
join --nocheck-order -j1 -t $'\t' -o 1.1,2.3   <(sort -t$'\t' -k1 mitogenomes_gb_accession) <(sort -t$'\t' -k1 ${nucl_gb_accession2taxid}) > genbank_taxid

# makeblastdb with cleaned_complete_partial_mitogenomes.fa and genbank_taxid meta file
makeblastdb -in cleaned_complete_partial_mitogenomes.fa -parse_seqids -blastdb_version 5 -taxid_map genbank_taxid -dbtype nucl -title complete_partial_mitogenomes.fa -out ${blast_dir}/complete_partial_mitogenomes/


# 4. BLAST of ASV's from ASV table against the MitoFish db with blastn ######################################################################################################################################

# set the percent id parameter for blast
perc_id=97

# output files from dada2 are in: /flash/RavasiU/Ayse/eDNA_workspace/dada2_out
# the files used for downstream steps are:
## asv_map
## asv_table.fasta
## asv_table_mapped.tab
## filtering_results.csv
## seqtab_nochim.csv

# mkdir and cd to the directory for blast work
query_asv_table=/flash/RavasiU/Ayse/eDNA_workspace/dada2_out/asv_table.fasta
path2db=${blast_dir}/complete_partial_mitogenomes

cd ${blast_dir%}
mkdir -p ${blast_dir}/${perc_id}

# load modules
ml bioinfo-ugrp-modules
ml ncbi-blast/2.10.0+

blastn -task blastn -db ${path2db} -query ${query_asv_table} -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" -out ${blast_dir}/${perc_id}/blast_perc_id${perc_id}.tab -perc_identity ${perc_id} -evalue 1e-2

