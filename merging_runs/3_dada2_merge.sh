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
   echo "Usage: $0 -a seqtab1 -b seqtab2 -c seqtab3 -o dada2_dir "
   echo -e "\t-a input file directory of seqtab 1"
   echo -e "\t-b input file directory of seqtab 2"
   echo -e "\t-c input file directory of seqtab 3"
   echo -e "\t-o output file directory for DADA2 script - dereplication of individual runs"
   exit 1 # Exit script after printing help
}

while getopts "a:b:c:o" opt
do
   case "$opt" in
      a ) seqtab_1="$OPTARG" ;;
      b ) seqtab_2="$OPTARG" ;;
      c ) seqtab_3="$OPTARG" ;;
      o ) dada2_dir="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# 2. Running DADA2 through the R scripts dada2.r #####################################################################################################################################################

# exporting paths to R and package libraries.
export PATH="/apps/unit/RavasiU/ayse/R-4.2.1/bin:$PATH"
export PATH="/home/a/ayse-oshima/R/x86_64-redhat-linux-gnu-library/4.1:$PATH"

/apps/unit/RavasiU/ayse/R-4.2.1/bin/Rscript --vanilla dada2_2_merge_seqtabs.R ${seqtab_1} ${seqtab_2} ${seqtab_3} ${dada2_dir}

## this script will pass arguments from command line to the  R script:

## directory for input seqtab # 1
#    path2seqtab_1 <- args[1]
## directory for input seqtab # 2
#    path2seqtab_2 <- args[2]
## directory for input seqtab # 3
#    path2seqtab_3 <- args[3]
## directory for this script's dada2 output files
#    dada2_dir <- args[4]

