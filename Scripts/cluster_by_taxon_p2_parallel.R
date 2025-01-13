
# Usage: Rscript cluster_by_taxon_p2.R [primer_pair] [infilename] [singleton_cutoff] [taxon_level] [rarefaction_level] [threads] [kingdom_column] [fasta_folder] [resume]

#			Description: Uses output from p2 to aggregate sequences at the desired taxonomic and rarefaction levels
#			Options:
#				1 [primer_pair] 		- the primer pair name (current default 'ITS1FLR3')
#				2 [infilename]			- .RData file containing phyloseq object from quick_and_dirty and make_phyloseq
#											Default: paste("minimap_",primer_pair,"/Phyloseq_Outfile_MIMP_",primer_pair,".RData",sep="")
#				3 [singleton_cutoff] 	- minimum number of sequences to carry forward from a given taxonomic assignment
#				4 [taxon_level]			- the taxonomic level at which to aggregate and subset sequences 
#				5 [rarefaction_level]	- the number of sequences to randomly subset from each taxonid
#				6 [threads]				- the number of threads to run
#				7 [kingdom_column]		- informs the script which column of the in file tax table contains kingdom information, for compatibility with samtools/minimap2 output from various databases
#				8 [fasta_folder]		- folder to which to save the taxon-aggregated, rarefied fasta files
#				9 [resume]				- a flag to tell the script to resume [?RESUME WHAT?] where it left off	
#			Output: cluster_by_taxon/fasta_folder/cluster_xxx.fa

# currently no support for taxon level other than six

args = commandArgs(trailingOnly=TRUE)
#print(args)

loadpackages <- function(){
	library('dplyr', quietly=T, verbose=F, warn.conflicts=F)
	library('Biostrings', quietly=T, verbose=F, warn.conflicts=F)
	library('phyloseq', quietly=T, verbose=F, warn.conflicts=F)
	library('doBy', quietly=T, verbose=F, warn.conflicts=F)
	library('foreach', quietly=T, verbose=F, warn.conflicts=F)
	library('doSNOW', quietly=T, verbose=F, warn.conflicts=F)
}

print("Loading R packages...")
suppressPackageStartupMessages(loadpackages())
print("R packages loaded successfully!")

filterdna <- function(x, strsearch) {
	x[which(sub(" rc","",names(x)) %in% strsearch)]
}

print("Reading in sequence and tax ids...")

output_table <- read.table("cluster_by_taxon/taxid_seqid.tsv", header=T)

tax_string<- read.table("cluster_by_taxon/taxids.tsv", header=T)

primer_pair <- args[1]
infilename <- args[2]
singleton_cutoff  <- as.integer(args[3])

taxon_level <- args[4]

rarefaction_level <- as.integer(args[5])

threads <- args[6]

kingdom_column <- as.numeric(args[7])

fasta_folder <- args[8]

resume <- args[9]

#primer_pair <- 'ITS54';infilename <- "minimap_ITS54/Phyloseq_Outfile_MIMP_ITS54.RData";singleton_cutoff  <- 20;taxon_level <- 'orderonly';rarefaction_level <- 1000;threads <- 96;kingdom_column <- 3;fasta_folder <- "taxon_cluster_orders";resume <- "F"

print(paste("primer_pair =", primer_pair))
print(paste("infilename  =", infilename ))
print(paste("singleton_cutoff =", singleton_cutoff))
print(paste("taxon_level =", taxon_level))
print(paste("rarefaction_level =", rarefaction_level))
print(paste("threads =", threads))
print(paste("kingdom_column =", kingdom_column))
print(paste("fasta_folder =", fasta_folder))
print(paste("resume =", resume))

print("Loading phyloseq object...")

load(infilename)

#load("minimap_2023-10-23-164538/Kelsey_phyloseq_workspace_ITS1FLR3.RData")

print("Filtering phyloseq object...")

if (kingdom_column != 1) {
	tax_table(psobj) <- tax_table(psobj)[,-c(1:(kingdom_column-1))]
	colnames(tax_table(psobj)) <- c('ta1','ta2','ta3','ta4','ta5','ta6','ta7')
}

print("Summarizing phyloseq object...")

abunds<- psobj %>%
  (function (x) data.frame(tax_table(x), numseqs=otu_table(x) %>% rowSums )) %>%
  summaryBy(numseqs ~ ta1 + ta2 + ta3 + ta4 + ta5 +ta6, data=., FUN=sum, var.names="abund") %>%
  left_join(tax_string) %>% na.omit

#primer_pair <- "ITS1FLR3"
#singleton_cutoff  <- 20
#taxon_level_cutoff <- 5
#rarefaction_level <- 1000
#max_taxon_grouping_size <- rarefaction_level*10

# unidentified sequences

###### THE FOLLOWING COMMENT MAY HAVE ALREADY BEEN DEPRECATED!

###########################################################
##  anything that's only identified to the family level  ##
##  will also be sub-clustered. anything not identified  ##
##  to order and above cannot be sub-clustered and will  ##
##  be treated separately/later...                       ##
###########################################################

########## DEFAULT::: GENUS OR FAMILY !!!!!!!!!!

## later on, will successively attempt to create alignments and clusters
## for things that aligned to accessions at higher taxon levels...

# things unidentified to family
# abunds %>% filter(ta5 == "unidentified")

if (taxon_level == "familygenus") {
	taxon_level_cutoff <- 5
	print(paste("Cutting off at", taxon_level, "(ta",taxon_level_cutoff,")..."))

	print(paste("Filtering things identified to taxon level ta", taxon_level_cutoff ,"..."))

	# things we don't want:

	# things not identified at least to family

	things_unidentified_to_taxonlevel <-
  	abunds[which(abunds[,paste("ta",taxon_level_cutoff,sep="")]=="unidentified"),'taxclusterid'] # things unidentified to family

} else if (taxon_level == "orderonly") {
	taxon_level_cutoff <- 4
	print(paste("Cutting off at", taxon_level, "(ta",taxon_level_cutoff,")..."))
	print(paste("Filtering things identified to taxon level ta", taxon_level_cutoff ,"..."))

	# things we don't want:
	# 1) things not identified at least to order
	# 2) things identified at least to family

	things_unidentified_to_taxonlevel <-
  	abunds[
		which(
			abunds[,paste("ta",taxon_level_cutoff,sep="")]=="unidentified" |  # 1) things not identified at least to order
			abunds[,paste("ta",taxon_level_cutoff+1,sep="")] != "unidentified"), # 2) things identified at least to family
		'taxclusterid']
} else if (taxon_level == "classonly") {
	taxon_level_cutoff <- 3
	print(paste("Cutting off at", taxon_level, "(ta",taxon_level_cutoff,")..."))
	print(paste("Filtering things identified to taxon level ta", taxon_level_cutoff ,"..."))

	# things we don't want:
	# 1) things not identified at least to class
	# 2) things identified at least to order

	things_unidentified_to_taxonlevel <-
  	abunds[
		which(
			abunds[,paste("ta",taxon_level_cutoff,sep="")]=="unidentified" |  # 1) things not identified at least to class
			abunds[,paste("ta",taxon_level_cutoff+1,sep="")] != "unidentified"), # 2) things identified at least to order
		'taxclusterid']
} else if (taxon_level == "phylumonly") {
	taxon_level_cutoff <- 2
	print(paste("Cutting off at", taxon_level, "(ta",taxon_level_cutoff,")..."))
	print(paste("Filtering things identified to taxon level ta", taxon_level_cutoff ,"..."))

	# things we don't want:
	# 1) things not identified at least to phylum
	# 2) things identified at least to class

	things_unidentified_to_taxonlevel <-
  	abunds[
		which(
			abunds[,paste("ta",taxon_level_cutoff,sep="")]=="unidentified" |  # 1) things not identified at least to phylum
			abunds[,paste("ta",taxon_level_cutoff+1,sep="")] != "unidentified"), # 2) things identified at least to class
		'taxclusterid']
} else if (taxon_level == "kingdomonly") {
	taxon_level_cutoff <- 1
	print(paste("Cutting off at", taxon_level, "(ta",taxon_level_cutoff,")..."))
	print(paste("Filtering things identified to taxon level ta", taxon_level_cutoff ,"..."))

	# things we don't want:
	# 1) things not identified at least to kingdom
	# 2) things identified at least to phylum

	things_unidentified_to_taxonlevel <-
  	abunds[
		which(
			abunds[,paste("ta",taxon_level_cutoff,sep="")]=="unidentified" |  # 1) things not identified at least to kingdom
			abunds[,paste("ta",taxon_level_cutoff+1,sep="")] != "unidentified"), # 2) things identified at least to phylum
		'taxclusterid']
} else {
	print(paste(taxon_level,"not supported!"))
	quit()
}


if (!(fasta_folder %in% dir('cluster_by_taxon'))) system(paste("mkdir cluster_by_taxon/",fasta_folder,sep=""))

tocyclethrough <- setdiff(tax_string$taxclusterid,things_unidentified_to_taxonlevel)

totaltax <- length(tax_string$taxclusterid)
numtax <- length(tocyclethrough)
numbar <- length(unique(output_table$barcode))


print(paste("Gathering and writing",numtax, "of",totaltax, "taxa from", numbar  ,"barcodes on",threads, "clusters to separate FASTA files..."))

cl <- makeCluster(spec=rep('localhost',
	#40),
	min(threads,numtax)),
	type="SOCK",
	outfile=paste("aggregate_fasta_log_",format(Sys.time(), "-%a%b%d-%H%M%S"),".csv", sep=""))

registerDoSNOW(cl)

pb <- txtProgressBar(min=0, max=numtax + 1, style=3)
#pb <- tkProgressBar(max=numtax)
progress <- function(n) setTxtProgressBar(pb, n)
#progress <- function() cat('.')
opts <- list(progress=progress)

x <- foreach (a=1:numtax,
#x <- foreach (a=1:40,
	#.combine='rbind',
	.final = function(x) NULL,
	.packages=c('dplyr', 'Biostrings'), .options.snow=opts) %dopar% {

#	general.error <- T
#	while(general.error) {tryCatch({
 
#for (a in 1:numtax) {
	all_seqs <- output_table %>% filter(taxclusterid == tocyclethrough[a])
	num_seqs <- dim(all_seqs)[1]
	if (
		num_seqs >= singleton_cutoff &
		!(
			as.logical(resume) & 
			(
				paste('cluster_', tocyclethrough[a], ".fa", sep="") %in%
				dir(
					paste('cluster_by_taxon/',fasta_folder,sep="")
					))
			)) {
		print(paste("Sampling taxid", tocyclethrough[a],paste0(abunds[abunds$taxclusterid==tocyclethrough[a],1:6],collapse=": "),"..."))
		if (num_seqs > rarefaction_level) {
			print(paste(rarefaction_level, "of", num_seqs, "sequences..."))

#			sample.error <- T
#			while(sample.error) {
#				sample.error <- tryCatch(
#					{
					samp_seqs <- all_seqs[sample(1:num_seqs, rarefaction_level),]#;F},
#					error = function(e) {print(paste("Sample error: taxid", tocyclethrough[a],paste0(abunds[abunds$taxclusterid==tocyclethrough[a],1:6],collapse=": ")));T},
#					warning = function(w) {F})
#			}
		} else {
			print(paste("all", num_seqs, "sequences..."))
			samp_seqs <- all_seqs
		}
		dna <- DNAStringSet()

		for (b in unique(samp_seqs$barcode)) {
			if ("all_filt_reorient.fasta" %in% dir(paste(b,primer_pair, sep="/"))) {
#				print(b)
				connection.error <- T
				while(connection.error) {
					connection.error <-
					tryCatch(
						{
						dna <- c(dna, filterdna(
							x=readDNAStringSet(paste(b,primer_pair,"all_filt_reorient.fasta", sep="/")),
							strsearch=samp_seqs$seqid))
						F},	error = function(e) {T}, warning = function(w) {F})
				}
			}

			# now write the dna for each taxonomic cluster and do a mafft alignment in another script!
			writeXStringSet(dna, paste('cluster_by_taxon/',fasta_folder,'/cluster_', tocyclethrough[a], ".fa", sep=""))

			#tocyclethrough[a]

			#c("Wrote", tocyclethrough[a], abunds[abunds$taxclusterid==tocyclethrough[a],1:6], a)
		}
	} #else {
	
		#tocyclethrough[a]
		
		#c("Skipped", tocyclethrough[a], abunds[abunds$taxclusterid==tocyclethrough[a],1:6], a)
	#}
	
#	F},	function(e) {print(paste("Error: taxid", tocyclethrough[a],paste0(abunds[abunds$taxclusterid==tocyclethrough[a],1:6],collapse=": ")));T}, warning = function(w) {F})}
}
close(pb)
stopCluster(cl)
#write.csv(x, paste("aggregate_fasta_log_",Sys.time(),".csv", sep=""),append=T)
