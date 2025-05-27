# this script parses the samtools output to make it readable in R, particularly in make_phyloseq.R

args = commandArgs(trailingOnly=TRUE)

# needs to change - quality score for alignment (mapq) is different than the quality score for the filtering...
qcut <- args[1]

minlen <- args[2]

##
#####

# Read SAM lines excluding headers
#print("Reading", getwd())
#getwd() %>% strsplit
sam_lines <- readLines("all_filt.samview.tsv")
sam_data_lines <- sam_lines[!grepl("^@", sam_lines)]

ntabs <- unlist(lapply(sam_data_lines, stringr::str_count, pattern = '\t'))
#sam_data_lines[order(ntabs, decreasing=T)][1:3]

# Convert to data.frame
#library(data.table)
#sam_out <- data.table::fread(text = sam_data_lines, header=FALSE, fill=T, sep='\t')

sam_out <- data.table::fread(text = sam_data_lines[order(ntabs, decreasing=T)], header=FALSE, fill=T, sep='\t')


# Add column names from SAM spec
colnames(sam_out)[1:11] <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR",
                      "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL")

sam_out <- sam_out[,c(1,2,3,4,5,10)]

#####

# hist(sam_out$MAPQ)
# length(unique(sam_out$QNAME))

#sam_out2 <- sam_out[-which(sam_out$mapq < qcut | nchar(sam_out$seq) < minlen),]
sam_out$seqlen <- nchar(sam_out$SEQ)
if (length(which(sam_out$seqlen < sam_out$minlen))>0) sam_out2 <- sam_out[-which(sam_out$seqlen < sam_out$minlen),]
sam_out2 <- sam_out2[which(sam_out2$FLAG %in% c(0,4,16,32,64)),]

# need to filter out higher scoring alignment per readid
sam_out3 <-
  dplyr::distinct(dplyr::select(dplyr::filter(dplyr::filter(dplyr::ungroup(dplyr::mutate(dplyr::group_by(sam_out2, QNAME),
                                             maxscore = max(MAPQ),
                                             maxlen   = max(seqlen))),
                              MAPQ == maxscore),seqlen == maxlen), -SEQ))

#relabund <- table(sam_out$tax_short)

#relabund2 <- sort(table(sam_out2$tax_short), T)
relabund3 <- sort(table(sam_out3$RNAME), T)

#check
# sam_out3 %>% with(table(readid)) %>% sort(decreasing = T) %>% head
# sam_out2 %>% filter(readid == '041dd021')
# sam_out %>% filter(readid == '041dd021')
# sam_out2 %>% filter(readid == 'cf7d10a2')
# sam_out %>% filter(readid == 'cf7d10a2')

# sam_out3 %>% filter(mapq < 5)
# length(unique(sam_out$readid))
# length(unique(sam_out3$readid))

# reads that did not map at cutoff

#setdiff(unique(sam_out$QNAME), unique(sam_out3$QNAME))

if(head(sort(table(sam_out3$QNAME), decreasing =T), 1)>1)
  {write.csv(data.frame(error="error! more than one unique maximum mapping per flag!"), "error.csv")
}else{
    
  write.csv(relabund3, paste("relabund_mapq", qcut ,".csv", sep=""), row.names=F)
  
  }

