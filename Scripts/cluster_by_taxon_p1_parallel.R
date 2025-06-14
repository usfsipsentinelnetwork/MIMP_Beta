# NOTE: currently only supports minimap2 output using UNITE 5-22-25

# Usage: Rscript cluster_by_taxon_p1.R [primer_pair] [infilename] [kingdom_column] [threads]
# required

#			Description: Tabulates seqids by taxon and outputs them to a file
			
#			Options:
#				1 [primer_pair] - the primer pair name (current default 'ITS1FLR3')
#				2 [infilename]	- .RData file containing phyloseq object from quick_and_dirty and make_phyloseq
#				3 [kingdom_column] - informs the script which column of the in file tax table contains kingdom information, for compatibility with samtools output from various databases (default
#				4 [threads]  		- number of threads to run
			
#			Output: cluster_by_taxon/taxid_seqid.tsv


#setwd("minimap_2023-10-23-164538")
#x <- read.csv("relabund_phred_q10.barcode01.csv")


print("Loading R packages...")

loadpackages <- function(){
	library('phyloseq', quietly=T, verbose=F, warn.conflicts=F)
	library('dplyr', quietly=T, verbose=F, warn.conflicts=F)
	library('progress', quietly=T, verbose=F, warn.conflicts=F)
	library('Biostrings', quietly=T, verbose=F, warn.conflicts=F)
	library('foreach', quietly=T, verbose=F, warn.conflicts=F)
	library('doSNOW', quietly=T, verbose=F, warn.conflicts=F)
}
suppressPackageStartupMessages(loadpackages())

print("R packages loaded successfully!")

args = commandArgs(trailingOnly=TRUE)

primer_pair <- args[1]

if (is.na(primer_pair))
	primer_pair <- "ITS1FLR3"

infilename <- args[2]

if (is.na(infilename))
	infilename <- paste("minimap_",primer_pair,"/Phyloseq_Outfile_MIMP_",primer_pair,".RData",sep="")
	
kingdom_column <- as.numeric(args[3])
threads <- as.numeric(args[4])

## infilename <- "minimap_2023-10-23-164538/Kelsey_phyloseq_workspace_ITS1FLR3.RData"

#print(primer_pair)
#print(infilename)
#quit()

print(paste("Loading", infilename))

load(infilename)
#psobj %>% tax_table %>% head

print(paste(infilename, "successfully loaded"))


tax_table(psobj) <- tax_table(psobj)[,-c(1:(kingdom_column-1))]
#colnames(tax_table(psobj)) <- c('ta1','ta2','ta3','ta4','ta5','ta6','ta7')
colnames(tax_table(psobj)) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

#print("Renamed columns")

# BOTH SHOULD WORK

### OLD VERSION, KEEP USING ta1-7
#genera_families <- tax_table(psobj)[,dim(tax_table(psobj))[2]-2:1] %>% as.data.frame %>% na.omit %>% distinct
#genera <- tax_table(psobj)[,dim(tax_table(psobj))[2]-1]  %>% as.data.frame %>% na.omit %>% distinct

### 5-22-25 FOR ALREADY-FORMATTED
genera_families <- tax_table(psobj)[,c("Family","Genus")] %>% as.data.frame %>% na.omit %>% distinct
genera <- tax_table(psobj)[,"Genus"]  %>% as.data.frame %>% na.omit %>% distinct

tax_string <- tax_table(psobj)[,-dim(tax_table(psobj))[2]] %>% as.data.frame %>% na.omit %>% distinct

# now make table with columns family, genus, barcode, sequence number
# and split the table up by genus or family and then cluster

barcodes <- grep("^barcode[0-9]+$", dir(), value=T, perl=T)

tax_string$taxclusterid = 1:(dim(tax_string)[1])

if (!('cluster_by_taxon' %in% dir())) system("mkdir cluster_by_taxon")

write.table(tax_string, "cluster_by_taxon/taxids.tsv", quote=F, row.names=F)

print(paste("Getting taxid-readnumber-barcode information..."))

niters <- length(barcodes)

chunks <- niters %/% threads

for (chunk in 0:chunks) {
#for (chunk in 1:chunks) {
    
	start0 <- 1 + threads*chunk
	stop0 <- min(start0 + threads - 1, niters)

	print(paste("Gathering sequences for", threads, "of",niters,"barcodes...(", round(100*chunk/(chunks+1)), "% finished )..."))
	cl <- makeCluster(spec=rep('localhost',threads), type="SOCK")
	registerDoSNOW(cl)
	
	# make a progress bar
	pb <- txtProgressBar(min=0, max=stop0-start0+1, style=3)
	progress <- function(n) setTxtProgressBar(pb, n)
	opts <- list(progress=progress)

	# cycle through the barcodes and taxa in parallel
	output_table_chunk <- foreach(a=start0:stop0, .combine="rbind", .options.snow=opts, .packages = c('dplyr','tidyr')) %dopar% {
		b <- barcodes[a]
		if (b %in% dir() & primer_pair %in% dir(b) & "all_filt.samview.tsv" %in% dir(paste(b,primer_pair,sep="/"))) {

			filename <- paste(b, "/", primer_pair, "/all_filt.samview.tsv", sep="")
			#x <- read.csv(filename, sep="\t", header=F, fill=T)

			#### May 22 2025 corrected ####
  			sam_lines <- readLines(filename)
  			sam_data_lines <- sam_lines[!grepl("^@", sam_lines)]
  			
  			ntabs <- unlist(lapply(sam_data_lines, stringr::str_count, pattern = '\t'))
  			#sam_data_lines[order(ntabs, decreasing=T)][1:3]
  			
  			# Convert to data.frame
  			#library(data.table)
  			#sam_out <- data.table::fread(text = sam_data_lines, header=FALSE, fill=T, sep='\t')
  			
  			x <- data.table::fread(text = sam_data_lines[order(ntabs, decreasing=T)], header=FALSE, fill=T, sep='\t')
  			
  			
  			# Add column names from SAM spec
  			colnames(x)[1:11] <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR",
  			                             "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL")
  			
  			x <- x[,c(1,2,3,4,5,10)]
  			x$seqlen <- nchar(x$SEQ)
  			if (length(which(x$seqlen < x$minlen))>0) {x <- x[-which(x$seqlen < x$minlen),]}
  			x <- x[which(x$FLAG %in% c(0,4,16,32,64)),]
  			
  			# need to filter out higher scoring alignment per readid
  			x <- x %>%
  			  group_by(QNAME) %>%
  			  mutate(maxscore = max(MAPQ), maxlen   = max(seqlen)) %>%
  			  ungroup %>%
  			  filter(MAPQ == maxscore) %>%
  			  filter(seqlen == maxlen) %>%
  			  select(-SEQ) %>%
  			  distinct
  			  
  			x <- x %>%
  			  separate(RNAME, remove = F, sep = '\\|', into = c("tax_unique","db","accession","accession_type","tax_string_full")) %>%
  			  mutate(tax_string =
  			           gsub("^(k__[A-Za-z_]+;p__[A-Za-z_]+;c__[A-Za-z_]+;o__[A-Za-z_]+;f__[A-Za-z_]+;g__[A-Za-z_]+)(;s__[A-Za-z_]+)$","\\1",tax_string_full, perl=T))
  			
			
			print(paste("Reading taxids for", b, "(", dim(x)[1], "sequence reads )")) # this doesn't do anything - need to add to the "opts" function or something

			# MODIFIED 5-22-25

			tax_string %>%
			  mutate(tax_string=paste("k__", Kingdom,
			                          ";p__", Phylum,
			                          ";c__", Class,
			                          ";o__", Order,
			                          ";f__", Family,
			                          ";g__", Genus, sep="")) %>%
			  left_join(x,.) %>% mutate(barcode=b) %>%select(taxclusterid, barcode, QNAME)
		}
	}
	
	close(pb)
	stopCluster(cl)
	print(paste("Writing taxid and barcode for", dim(output_table_chunk)[1], "reads to file..."))

	if (chunk == 0) {
	  write.table(output_table_chunk, "cluster_by_taxon/taxid_seqid.tsv", quote=F, row.names=F)
	} else {
	  write.table(output_table_chunk, "cluster_by_taxon/taxid_seqid.tsv", quote=F, row.names=F, append=T, col.names=F)
	}
}

final_output <- read.table("cluster_by_taxon/taxid_seqid.tsv", header=1) %>% distinct
write.table(final_output, "cluster_by_taxon/taxid_seqid.tsv", quote=F, row.names=F)
#with(read.table("taxid_seqid.tsv", header=1),write.table(cbind(seqid,seqid), 'seqids.names', quote=F, row.names=F, col.names=F))