# need development version of dplyr ## not any more

# Usage: Rscript make_phyloseq.R [primer_pair*] [db*] [outfilename] [silva_path]
# * required

#setwd("C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Projects/MinIon/KelseyTobin_results")

#	Usage: Rscript make_phyloseq.R [primer_pair*] [db*] [outfilename] [silva_path]

#	Description: Uses preliminary quick and dirty minimap2 to create an OTU and taxon table to make a phyloseq object

#	Options (* required)
#		1 [primer_pair] - the primer pair name
#		2 [db]			- the database (currently supports =UNITE and =SILVA)
#		3 [outfilename] - name of .RData output file with phyloseq object
#		4 [silva_path]  - ***required if db=SILVA

#	Output: outfilename.RData

args = commandArgs(trailingOnly=TRUE)

primer_pair <- args[1] #

db <- args[2] # name of database ('UNITE' or 'SILVA')

# options
# UNITE
# SILVA - if SILVA must supply path

# this isn't working
outfilename <- args[3]

if (is.na(outfilename))
	outfilename <- paste("Phyloseq_Outfile_MIMP_",primer_pair,".RData",sep="")

if (db == "SILVA") silva_path <- args[4] #"C:/Users/GeoffreyWilliams/OneDrive - USDA/Project_folders_Personal/Projects/MinIon/MIMP/ref_dbs/taxmap_slv_ssu_ref_nr_138.1.txt" #args[3]

setwd(paste("minimap", primer_pair, sep="_"))

library(dplyr)    ### NO LONGER NEEDED: development version has join_by
library(tidyr)
library(phyloseq)

dd <- dir()

output_matrix <- NULL

for (filename in grep("^relabund.+\\.csv",dd,value=T,perl=T)) {
  
  current_table <- read.csv(filename)
  current_table$barcode <- strsplit(filename, split="\\.")[[1]][2]
  
  output_matrix <- rbind(output_matrix, current_table)
}

if (db == "SILVA") {
  
  taxdb <- read.csv(silva_path, sep="\t")
  
  output_matrix <- separate(output_matrix, Var1, into=c("primaryAccession", "start", "stop"), sep="\\.")

  #output_matrix$start <- as.numeric(output_matrix$start)
  #output_matrix$stop <- as.numeric(output_matrix$stop)
    
  output_matrix <- left_join(output_matrix, taxdb, by=c("primaryAccession"))
  output_matrix_cast <- reshape2::acast(output_matrix, primaryAccession + path + organism_name ~ barcode, value.var='Freq')
  output_matrix_cast[which(is.na(output_matrix_cast))]<-0
  
  df1<- data.frame(taxon_names = rownames(output_matrix_cast))
  row_names_long <-
    df1 %>%
    separate(
      taxon_names,
      into=c("Accession1","Taxonomy","Species"), sep="_") %>%
    separate(Taxonomy,
             into=c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";")
  
  
} else if (db == "UNITE") {
  output_matrix_cast <- reshape2::acast(output_matrix, Var1 ~ barcode, value.var='Freq')
  output_matrix_cast[which(is.na(output_matrix_cast))]<-0
  
  df1<- data.frame(taxon_names = rownames(output_matrix_cast))
  row_names_long <-
    df1 %>%
    separate(
      taxon_names,
      into=c("Species","Accession1","Accession2","type_of_seq","Taxonomy"), sep="\\|") %>%
    mutate(Taxonomy=gsub("[kpcofgs]__","",Taxonomy)) %>%
    separate(Taxonomy,
             into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=";")  
}

row.names(output_matrix_cast) <- 1:dim(output_matrix_cast)[1]

row.names(row_names_long)<-1:dim(row_names_long)[1]

OTU <- otu_table(output_matrix_cast, taxa_are_rows=T)
TAX <- tax_table(row_names_long[,-1])

#taxa_names(OTU)

taxa_names(TAX) <- taxa_names(OTU)
psobj<-phyloseq(OTU, TAX)

save(psobj, file=outfilename)
