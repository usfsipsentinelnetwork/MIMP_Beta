#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("R.utils")
library("dplyr")

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  args[1] = "all_filt.denovo.fasta"
}

#args[1]
#l <- countLines("all_filt.denovo.fasta")/2
#l <- countLines(args[1])/2
l <- min(500,countLines(args[1])/2)

#args[2]
#align<-read.table("all_filt_denovo.maf.cons.data", row.names=NULL)
align<-read.table(args[2], row.names=NULL)

bases<-names(align)[2:5]

align$max_score<-sapply(rownames(align), FUN=function(x) max(align[x,bases]))
align$proportion <- align$max_score/l
align$over50 <- align$proportion > .50

# second highest score:
### innermost is index of 2nd highest score: order(align[x,c("a","c","g","t")], decreasing=T)[2]
### next innermost is the 2nd highest score, obtained as following:
###### index of second highest score is applied to scores, which are align[x,c("a","c","g","t")]
align$second_highest_score <- sapply(
	rownames(align),
	FUN=function(x) #as.numeric(align[x,bases])[as.numeric(order(align[x,bases])[3])])


as.numeric(
	align[rownames(align)[1],bases]
	)[
		as.numeric(
			order(
				unlist(align[rownames(align)[1],bases])
			)[3])]
)

# now calculate ratio of highest and second highest, and identify bases where it is greater than 2 (double)
align$proportion_over_2 <- align$max_score/align$second_highest_score >=1.5

#align[x,c("a","g","c","t")][]
# now base of highest score is obtained in similar fashion
### innermost is index of highest score: order(align[x,c("a","g","c","t")], decreasing=T)[1]
### next innermost is the base with highest score, obtained as following:
### index of highest score is applied to the bases, which are c("a","g","c","t")
align$highest_scoring_base <- sapply(
	rownames(align),
	FUN=function(x) bases[as.numeric(order(
							unlist(align[x,bases])
						)[4])])

# now find the indices where above is the case, but base was called as an X (for some reason?)
bases_to_replace <- which(align$row.names=="X" & align$proportion_over_2)

# now replace those bases with
align$row.names[bases_to_replace] <- align$highest_scoring_base[bases_to_replace]

cons<-paste0(align[which(align$over50), "row.names"],collapse="")
score1<- with(align[which(align$over50),], sum(proportion)/length(proportion))

seqname <- strsplit(getwd(),"/")[[1]] %>% (function(x) paste(x[length(x)-1], x[length(x)], args[2]))

conn1<-file(paste(as.character(gsub(" ",".",seqname)),"cons.denovo.python.fasta",sep="."))
writeLines(c(paste(">", seqname, " score=", score1, " reads=", l, sep=""), as.character(gsub("X","",cons))), conn1)

close(conn1)

#conn2<-file(paste(seqname,"cons.python.score",sep="."))
#writeLines(as.character(score1), conn2)
#close(conn2)

write.csv(align[which(align$over50),], paste(as.character(gsub(" ",".",seqname)),"cons.denovo.python.data",sep="."), row.names=F)
