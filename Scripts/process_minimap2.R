# this script parses the samtools output to make it readable in R, particularly in make_phyloseq.R

args = commandArgs(trailingOnly=TRUE)

qcut <- args[1]

minlen <- args[2]

##
sam_out <- read.csv("all_filt.samview.tsv", sep="\t", header=F, fill=T)[,c(1,2,3,5,10)]

names(sam_out) <- c("readid", "flag", "taxid", "phredscore", "seq")

#sam_out$tax_short <- sapply(sam_out$taxid, function (x) ifelse(grepl(" ",x),paste(paste(strsplit(x," ")[[1]][2:length(strsplit(x, " ")[[1]])],collapse="_"),strsplit(x," ")[[1]][1],sep="|"),paste(strsplit(x, "\\|")[[1]][1],strsplit(x, "\\|")[[1]][2],strsplit(x, "\\|")[[1]][3],sep="|")))

sam_out2 <- sam_out[-which(sam_out$phredscore < qcut & nchar(sam_out$seq) < minlen),]

sam_out2 <- sam_out2[which(sam_out2$flag %in% c(0,4,16,32,64)),]

#relabund <- table(sam_out$tax_short)

#relabund2 <- sort(table(sam_out2$tax_short), T)
relabund3 <- sort(table(sam_out2$taxid), T)

#write.csv(relabund2, paste("relabund_phred_q", qcut ,".csv", sep=""), row.names=F)
write.csv(relabund3, paste("relabund_phred_q", qcut ,".csv", sep=""), row.names=F)

