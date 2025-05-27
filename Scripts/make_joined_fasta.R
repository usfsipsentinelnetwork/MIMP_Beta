
# usage:
#
# shell script should already have catted the barcode files together
#
# Rscript make_joined_fasta.R [primer_pair] [taxidfile] [taxidseqidfile] [aggregated_fasta_file]

# setwd("C:/Users/readi/OneDrive/Project_folders_Personal/USDA2021to24/Projects/MinIon/Argentina_Backup/MinIon2022_forLucia_Feb2024/LM1-recalled-demultiplexed-80")

print("Loading R packages...")

loadpackages <- function(){
  library('phyloseq', quietly=T, verbose=F, warn.conflicts=F)
  library('dplyr', quietly=T, verbose=F, warn.conflicts=F)
  #library('progress', quietly=T, verbose=F, warn.conflicts=F)
  library('Biostrings', quietly=T, verbose=F, warn.conflicts=F)
  library('dplyr')
  library('tidyr')
#  library('foreach', quietly=T, verbose=F, warn.conflicts=F)
#  library('doSNOW', quietly=T, verbose=F, warn.conflicts=F)
}
suppressPackageStartupMessages(loadpackages())

print("R packages loaded successfully!")

args = commandArgs(trailingOnly=TRUE)

primer_pair <- args[1]

if (is.na(primer_pair))
  primer_pair <- "ITS54"

infilename <- args[2]

if (is.na(taxidfile))
  taxidfile <- "cluster_by_taxon/taxids.tsv"

infilename <- args[3]

if (is.na(taxidfile))
  taxidseqidfile <- "cluster_by_taxon/taxid_seqid.tsv"

infilename <- args[4]

if (is.na(taxidfile))
  aggregated_fasta_file <- "all_sequences.fasta"

taxid <- read.table(taxidfile, sep = ' ', header=T)
taxidseqid <- read.table(taxidseqidfile, sep = ' ', header=T)
#head(taxid)
#head(taxidseqid)

dna <- readDNAStringSet(aggregated_fasta_file)
#dna %>% as.data.frame

#str(dna)

#nms <- data.frame(names_orig = names(dna)[1:5]) %>% separate(names_orig, into = c('seqid','rc'))

#alphabet(dna)
#dna[[1]]
#as.character(dna, use.names=F)

#dna_names_joined1to5 <-
#  data.frame(names_orig = names(dna)[1:5]) %>%
#  separate(names_orig, into = c('seqid','rc')) %>%
#  cbind(seq = as.character(dna, use.names=F)[1:5]) %>%
#  select(-rc) %>%
#  left_join(taxidseqid) %>%
#  left_join(taxid)
  
#table(taxidseqid$seqid) %>% sort(decreasing =T) %>% head

#filter(taxidseqid, seqid=="A_331_3637")seqid


dna_names_joined <-
  data.frame(names_orig = names(dna)) %>%
  separate(names_orig, into = c('seqid','rc'), sep=" ") %>%
  cbind(seq = as.character(dna, use.names=F)) %>%
  left_join(taxidseqid, by = join_by(seqid == QNAME)) %>% #dim
  left_join(taxid) %>%
  select(barcode, seqid, taxclusterid, Kingdom, Phylum, Class, Order, Family, Genus, seq)

write.csv(dna_names_joined, "joined_fasta_barcodes_clusters.csv", row.names = F)

dna_names_joined_seqobj <- dna

names(dna_names_joined_seqobj)<-with(dna_names_joined, paste(barcode, seqid, taxclusterid, Kingdom, Phylum, Class, Order, Family, Genus))

names(dna_names_joined_seqobj)

writeXStringSet(dna_names_joined_seqobj, "joined_fasta_barcodes_taxonomy.fasta")