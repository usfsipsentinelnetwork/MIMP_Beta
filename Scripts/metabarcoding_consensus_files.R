library("R.utils")

otu_table<- read.csv("outfile.an.list", sep="\t")

#if ( length(dir("all_filt_reorient.an.0.1.fasta"))==0 )
#{
#	cutoff <- otu_table$label[dim(otu_table)[1]]
#} else {
#	cutoff="0.1"
#}

cutoff<-0.51

filename <- paste("outfile.an.", cutoff, ".fasta", sep="")
filename2 <- paste("outfile.an.", cutoff, ".notabs.fasta", sep="")
filename3 <- paste("outfile.an.", cutoff, ".temp.fasta", sep="")
filename4 <- paste("outfile.an.", cutoff, ".temp2.fasta", sep="")

system(paste("sed 's/\t/_/g' ", filename, " > ", filename2, sep=""))

label_otus <- otu_table[otu_table$label==cutoff, 1:2]

x <- 0

# make separate files for each OTU

num_otus <- as.numeric(label_otus["numOtus"])

system(paste("cat", filename2, ">", filename3))
system(paste("cat", filename3, ">", filename4))

for (otu in 1:num_otus) {
	filename_otu <- paste(names(otu_table)[2+otu], ".fasta", sep="")
	
	# take first x sequences (from current otu) of filename4 and put in filename5
	# x is meant to be size of OTU
	x <- length(strsplit(otu_table[otu_table$label==cutoff, otu+2], ",")[[1]])
	
	if (x < 3) break
	
	system(paste("head -n ", x*2, " ", filename3, " > ", filename_otu, sep=""))
	
	# cut out first x sequences of filename3, temporarily store in filename4
	total_lines <- as.numeric(countLines(filename3))

	system(paste("tail -n ", total_lines - x*2, " ", filename3, " > ", filename4, sep=""))
	# now replace filename3 with filename4
	system(paste("cat", filename4, ">", filename3))
}