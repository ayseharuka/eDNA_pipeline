#!/usr/bin/env Rscript
# DADA2 for trimming off primers and merging the reads
.libPaths("/home/a/ayse-oshima/R/x86_64-redhat-linux-gnu-library/4.1")
args = commandArgs(trailingOnly=TRUE)

# directory for input cutadapt files
cutadapt_dir <- args[1]
# directory for dada2 output files
dada2_dir <- args[2]
# path to and name of the seqence table 
seqtab_dir_name <- args[3]

library(dplyr)
library(tidyr)
library(tibble)
library(dada2)
library(Biostrings)
library(ShortRead)
library(stats)

# input path: cutadapt processed fastq files
# path.input <- file.path("/flash/RavasiU/Ayse/eDNA_workspace/cutadapt_out")
path.input <- file.path(cutadapt_dir)
path.input
print("the input path for dada2:\n")
print(countFastq(path.input))

# the output folder for dada2
# path.output <- file.path("/flash/RavasiU/Ayse/eDNA_workspace/dada2_out")
path.output <- file.path(dada2_dir)
if(!dir.exists(path.output)) dir.create(path.output)
print("the output path for dada2:\n")
print(path.output)

cutFs <- sort(list.files(path.input, pattern = "R1.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.input, pattern = "R2.fastq", full.names = TRUE))
print("printing the R1 files:\n")
print(cutFs)
print("printing the R2 files:\n")
print(cutRs)

print("learnErrors start")
errF <- learnErrors(cutFs, multithread=TRUE)
errR <- learnErrors(cutRs, multithread=TRUE)
print("learnErrors end")

# i could just save these as jpeg?
# plotErrors(errF, nominalQ = TRUE)
# plotErrors(errR, nominalQ = TRUE)

print("Dereplication start\n")
derepFs <- derepFastq(cutFs, verbose = TRUE)
derepRs <- derepFastq(cutRs, verbose = TRUE)
print("Dereplication - ok! \n")

# example: A1_trimmed_R1.fastq
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
# head(sample.names)
print("sample names:") 
sample.names

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct Sequence Table, ASV table!
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, seqtab_dir_name) 

