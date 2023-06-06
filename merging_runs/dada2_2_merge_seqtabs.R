#!/usr/bin/env Rscript

# this script is for merging sequence tables obtained previously from three separate runs

# DADA2 for trimming off primers and merging the reads
.libPaths("/home/a/ayse-oshima/R/x86_64-redhat-linux-gnu-library/4.1")

args = commandArgs(trailingOnly=TRUE)

# directory for input seqtab # 1
path2seqtab_1 <- args[1]
# directory for input seqtab # 2
path2seqtab_2 <- args[2]
# directory for input seqtab # 3
path2seqtab_3 <- args[3]

# directory for this script's dada2 output files
dada2_dir <- args[4]

library(dplyr)
library(tidyr)
library(tibble)
library(dada2)
library(Biostrings)
library(ShortRead)
library(stats)

# the output folder for dada2
path.output <- file.path(dada2_dir)
if(!dir.exists(path.output)) dir.create(path.output)
print("the output path for dada2:\n")
print(path.output)

# read in the seqtabs of multiple runs
seqtab1 <- readRDS(path2seqtab_1)
seqtab2 <- readRDS(path2seqtab_2)
seqtab3 <- readRDS(path2seqtab_3)
# merge the seqtabs into one
seqtab.all <- mergeSequenceTables(seqtab1, seqtab2, seqtab3)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) 
print("seqtab.nochim, the asv table:")
head(seqtab.nochim[, 5])

table(nchar(getSequences(seqtab.nochim)))

sum(seqtab.nochim)/sum(seqtab.all)

# write it if for the first time, read it in if not.
path.seqtab.nochim <- file.path(path.output, "seqtab_nochim.csv")
path.seqtab.nochim
write.csv(seqtab.nochim, path.seqtab.nochim, append = FALSE)
# seqtab.nochim <- read.csv("plots/seqtab_nochim.csv", header = TRUE, row.names = 1)

# transpose the seqtab.nochim as a dataframe so it is in the acceptable format for downstream
asv_table <- as.data.frame(t(as.matrix(seqtab.nochim))) %>%
  rownames_to_column("seq")
print("asv_table:")
head(asv_table)

# make ASV_map 
# take seq column and put into new table, add #ID of ASV1, ASV2, ASV3 ...
asv_seqs <- asv_table["seq"]
print("asv_seqs:")
head(asv_seqs)

asv_ids <- data.frame(paste0(rep("ASV", nrow(asv_seqs)), 1:nrow(asv_seqs)))
colnames(asv_ids) <- "#ID"
print("asv_ids:")
head(asv_ids)
asv_map <- cbind(asv_ids, asv_seqs)

print("asv_map:")
head(asv_map)

# saved both as csv and tab
path.asv.map.csv <- file.path(path.output, "asv_map.csv")
write.csv(asv_map, "path.asv.map")

path.asv.map.table <- file.path(path.output, "asv_map")
write.table(asv_map, path.asv.map.table, sep = "\t",
            quote = FALSE,
            row.names = FALSE, col.names = TRUE,
            append = FALSE)

asv_table_mapped <- inner_join(asv_map, asv_table, by = "seq") %>%
  select(-c("seq"))
print("asv_table_mapped:")
head(asv_table_mapped)

path.asv.table.mapped.tab <- file.path(path.output, "asv_table_mapped.tab")
path.asv.table.mapped.tab
# needs to be tab-delimited (with ext .tab) for the LCA script
write.table(asv_table_mapped, file = path.asv.table.mapped.tab, sep = "\t",
            quote = FALSE, # important!
            row.names = FALSE, col.names = TRUE, 
            append = FALSE)

# also need to create a fasta file of: 
# >ASV1
# ATGC.....
# >ASV2
# ATGC.....

# with BioStrings
path.asv.fasta <- file.path(path.output, "asv_table.fasta")
asv_fasta <- DNAStringSet(asv_map$seq)
names(asv_fasta) <- asv_map$`#ID`
writeXStringSet(asv_fasta, filepath = path.asv.fasta, format = "fasta")

#####################################################################################
# We now inspect the the number of reads that made it through each step in the pipeline to verify everything worked as expected.

#getN <- function(x) sum(getUniques(x))
#track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)

#colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
#rownames(track) <- sample.names
#print("tracking read counts through steps:")
#head(track)

#path.filter.results <- file.path(path.output, "filtering_results.csv")
#write.csv(track, path.filter.results)
#####





