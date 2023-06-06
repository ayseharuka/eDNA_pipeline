#!/bin/bash -l
#SBATCH --partition=compute
#SBATCH --mem=50G
#SBATCH --cpus-per-task=50
#SBATCH --mail-user=ayse-oshima@oist.jp
#SBATCH --mail-type=FAIL,END
#SBATCH --time=45

#### Synopsis:
# --> 1. Cutadapt trimming off primers, adapter seqs, filtering out N-containing reads, quality cut off etc.
#     2. Merging, dereplication, denoising to get ASV table with DADA2 (QIIME plug-in)
#     3. BLAST the ASV table against the makblastdb'd DB made from MitoFish seqs
#     4. LCA script

helpFunction()
{
   echo ""
   echo "Usage: $0 -i input_directory -o output_directory"
   echo -e "\t-i input file directory of the fastq files"
   echo -e "\t-o output file directory for CutAdapt outputs"
   exit 1 # Exit script after printing help
}

while getopts "i:o:" opt
do
   case "$opt" in
      i ) input_dir="$OPTARG" ;;
      o ) output_dir="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# show what were the arguments
echo ${input_dir}
echo ${output_dir}

# the path to the output directory, mkdir it under the working directory where the script is run.
mkdir -p ${output_dir}

# Load the necessary modules (loaded in the order of their usage):
ml bioinfo-ugrp-modules DebianMed
ml cutadapt/3.2-2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
##for array jobs, use $SLURM_ARRAY_TASK_ID, and submit as "sbatch --array=0-31, --array=1,4,10..."
tid=$SLURM_ARRAY_TASK_ID

# the delimiter I will use to awk for the basename
FOR=_S

# quality filtering threshold
qv=10

# create the cutadapt sanity check file with the header (id     r1_before       r1_after        r2_before       r2_after)
echo -e "id\tr1_before\tr1_after\tr2_before\tr2_after" > ${output_dir}/cutadapt_before_after.txt

#### Loop for grouping R1 and R2 files.

# gets sample base names (id): A1, A2, A3, A-Field etc.
for id in `ls ${input_dir} |awk -F ${FOR%} '{print $1}' |sort |uniq`
do

# grep for the base names (id) and R1 and R2 to get list of R1 and R2 fastq file names
R1=`ls ${input_dir} |grep $id |grep R1`
R2=`ls ${input_dir} |grep $id |grep R2`

# CutAdapt arguments:
# 1-4: processsing
## 1. -u length (+) from the beginning, (-) from the end
## 2. -q quality cut off
## 3. -a -g -b adapter/primer
## 4. -l length (+) from the end, (-) from the beginning (why not keep it consistent as the 1. step length -u? confusing.)
# 5-? filtering of processed reads: too short or untrimmed reads
## 5. -m discard reads shorter than this -M discard reads longer than this
## 6. --max-n discard reads with count of N (ambiguous base) greater than the given value. --> set to 0 to remove reads with N>0.

# Cutadapt command with arguments
cutadapt --cores ${OMP_NUM_THREADS} -u 6 -q ${qv},${qv} -g GTCGGTAAAACTCGTGCCAGC -g GCCGGTAAAACTCGTGCCAGC -g RGTTGGTAAATCTCGTGCCAGC -a CAAACTGGGATTAGATACCCCACTATG -a CAAACGGGGATTAGACACCCTCCTATG -a CAAACTAGGATTAGATACCCCACTATGC -G CATAGTGGGGTATCTAATCCCAGTTTG -G CATAGGAGGGTGTCTAATCCCCGTTTG -G GCATAGTGGGGTATCTAATCCTAGTTTG -A GCTGGCACGAGTTTTACCGAC -A GCTGGCACGAGTTTTACCGGC -A GCTGGCACGAGATTTACCAACY -n 2 -m 50 --max-n 0 -o ${output_dir}/${id}_trimmed_R1.fastq -p ${output_dir}/${id}_trimmed_R2.fastq ${input_dir}/${R1} ${input_dir}/${R2}

echo OK-Cutadapt

# append the the number of reads before and after the cutadapt to the sanity check file

# input fastq files:
# ${input_dir}/${R1}
# ${input_dir}/${R2}

# output fq files:
# ${output_dir}/${id}_trimmed_R1.fastq
# ${output_dir}/${id}_trimmed_R2.fastq

before_r1=`zgrep "@" ${input_dir}/${R1} | wc -l`
before_r2=`zgrep "@" ${input_dir}/${R2} | wc -l`
after_r1=`grep "@" ${output_dir}/${id}_trimmed_R1.fastq | wc -l`
after_r2=`grep "@" ${output_dir}/${id}_trimmed_R2.fastq | wc -l`
echo -e $id$'\t'$before_r1$'\t'$after_r1$'\t'$before_r2$'\t'$after_r2 >> ${output_dir}/cutadapt_before_after.txt

done

