#!/bin/bash -l
#SBATCH --partition=compute
#SBATCH --mem=50G
#SBATCH --cpus-per-task=50
#SBATCH --mail-user=ayse-oshima@oist.jp
#SBATCH --mail-type=FAIL,END
#SBATCH --time=45

#### Synopsis:
#     1. Cutadapt trimming off primers, adapter seqs, filtering out N-containing reads, quality cut off etc.
# --> 2. Merging, dereplication, denoising to get ASV table with DADA2 (QIIME plug-in)
#     3. BLAST the ASV table against the makblastdb'd DB made from MitoFish seqs
#     4. LCA script

helpFunction()
{
   echo ""
   echo "Usage: $0 -i input_directory -o output_directory"
   echo -e "\t-i input file directory of the CutAdapt outputs"
   echo -e "\t-o output file directory for DADA2 script - dereplication of individual runs"
   echo -e "\t-s path to and name of the seqtable from this individual run"
   exit 1 # Exit script after printing help
}

while getopts "i:o:s:" opt
do
   case "$opt" in
      i ) cutadapt_dir="$OPTARG" ;;
      o ) dada2_dir="$OPTARG" ;;
      s ) seqtab_dir_name="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# make the directory before calling the R script so the rds file can be saved there.
mkdir -p ${dada2_dir}

# 2. Running DADA2 through the R scripts dada2.r #####################################################################################################################################################

# exporting paths to R and package libraries.
export PATH="/apps/unit/RavasiU/ayse/R-4.2.1/bin:$PATH"
export PATH="/home/a/ayse-oshima/R/x86_64-redhat-linux-gnu-library/4.1:$PATH"

/apps/unit/RavasiU/ayse/R-4.2.1/bin/Rscript --vanilla dada2_1_seqtab.R ${cutadapt_dir} ${dada2_dir} ${seqtab_dir_name}

## this script will pass arguments from command line to the  R script:

## in the R script:

## 1. directory for input cutadapt files
#     cutadapt_dir <- args[1]

## 2. directory for dada2 output files
#     dada2_dir <- args[2]

## 3. path to and name of the seqence table
#     seqtab_dir_name <- args[3]

