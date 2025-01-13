#!/bin/bash
#	Usage: sh aggregate_minimap_bytaxon.sh [-p primer_pair] [-f phyloseq filename] [-l taxon level for clustering] [-h print this help message] [-s singleton_cutoff] [-r rarefaction_level] [-k kingdom_column] [-t threads] [-T tabulate only] [-A aggregate only] [-R resume] [-F fasta_folder]
#
#	Description: Uses output from preliminary quick and dirty minimap2 to do run the following R scripts
#		(supplying the options given to it to the R scripts); descriptions and specs for R scripts below
#		5a. cluster_by_taxon_p1_parallel.R
#			Output: cluster_by_taxon/taxid_seqid.tsv	
#		5b. cluster_by_taxon_p2_parallel.R	
#			Output: cluster_by_taxon/fasta_folder/cluster_xxx.fa
#	Options:
#		-p primer_pair='ITS54'		- the primer pair name (current default 'ITS1FLR3')
#		-f psfilename=''			- .RData file containing phyloseq object from quick_and_dirty and make_phyloseq
#		-l taxonlevel='genus'		- the taxonomic level at which to aggregate and subset sequences 
#		-s singleton_cutoff=20		- minimum number of sequences to carry forward from a given taxonomic assignment
#		-r rarefaction_level=1000	- the number of sequences to randomly subset from each taxonid
#		-k kingdom_column=3			- informs the script which column of the in file tax table contains kingdom information, for compatibility with samtools/minimap2 output from various databases
#		-t threads=8				- the number of threads to run
#		-T tabulate_only='F'		- whether to only run cluster_by_taxon_p1_parallel.R, which produces a seqid-by-taxid table cluster_by_taxon/taxid_seqid.tsv
#		-A aggregate_only='F'		- whether to only run cluster_by_taxon_p2_parallel.R, which does the aggregation by taxon and rarefaction
#		-R fasta_folder='taxon_cluster' - folder to which to save the taxon-aggregated, rarefied fasta files
#		-F resume='F'				- a flag to tell the script to resume [?RESUME WHAT?] where it left off	
#		-L session_log_file="[-p primer_pair].log" 
#				name of log file where you are keeping track of all the commands
#				you are running and with which parameters (automatically supplies
#				one based on primer pair)
#	Output (via R scripts):
#		cluster_by_taxon/taxid_seqid.tsv
#		cluster_by_taxon/fasta_folder/cluster_xxx.fa

#default parameter values
primer_pair='ITS54'
psfilename=''
taxonlevel='genus'
singleton_cutoff=20
rarefaction_level=1000
kingdom_column=3
threads=8
tabulate_only='F'
aggregate_only='F'
fasta_folder='taxon_cluster'
resume='F'
session_log_file=""

print_usage(){
	echo "Usage: sh aggregate_minimap_bytaxon.sh [-p primer_pair] [-f phyloseq filename] [-l taxon level for clustering] [-h print this help message] [-s singleton_cutoff] [-r rarefaction_level] [-k kingdom_column] [-t threads] [-T tabulate only] [-A aggregate only] [-R resume] [-F fasta_folder]"
}
while getopts 'p:f:l:s:r:k:t:TARF:L:h' OPTION; do
	case "$OPTION" in
		p)
			primer_pair="$OPTARG"
			;;
		f) psfilename="$OPTARG" ;;
		l) taxonlevel="$OPTARG" ;;
		s) singleton_cutoff="$OPTARG" ;;
		r) rarefaction_level="$OPTARG" ;;
		k) kingdom_column="$OPTARG" ;;
		t) threads="$OPTARG" ;;
		T) tabulate_only='T' ;;
		A) aggregate_only='T' ;;
		R) resume='T' ;;
		F) fasta_folder="$OPTARG" ;;
		L) session_log_file="$OPTARG" ;;
		h)
			print_usage
			exit ;;
		*)
		exit 1 ;;
	esac
done

current_time=$(date "+%Y.%m.%d-%H.%M.%S")
file_name="align-and-cluster-subtaxon-log"
file_name_stamped=$file_name.$current_time.txt

# check for ps object file name
if [[ $psfilename == "" ]]
then
	psfilename="minimap_$primer_pair/Phyloseq_Outfile_MIMP_$primer_pair.RData"
fi

# if session log file name not supplied, create one (e.g., ITS54.log)
if [[ $session_log_file == "" ]]; then
	log_file=$primer_pair.log
fi

# if session log file by that name does not yet exist, create a file
if ! [[ -f $session_log_file ]]; then
	>$session_log_file
fi

echo "running aggregate_minimap_bytaxon.sh -p $primer_pair -f $psfilename -l $taxonlevel -s $singleton_cutoff -r $rarefaction_level -k $kindgom_column -t $threads -T $tabulate_only -A $aggregate_only -F $fasta_folder -R $resume"
printf "running aggregate_minimap_bytaxon.sh -p $primer_pair -f $psfilename -l $taxonlevel -s $singleton_cutoff -r $rarefaction_level -k $kindgom_column -t $threads -T $tabulate_only -A $aggregate_only -F $fasta_folder -R $resume" > $file_name_stamped

# next output command and options to the session log file:
echo $current_time: >> $session_log_file
echo "sh aggregate_minimap_bytaxon.sh -p $primer_pair -f $psfilename -l $taxonlevel -s $singleton_cutoff -r $rarefaction_level -k $kindgom_column -t $threads -T $tabulate_only -A $aggregate_only -F $fasta_folder -R $resume" >> $session_log_file
echo "..." >> $session_log_file

echo "Using $primer_pair"
echo "Using phyloseq object file $psfilename"
#echo $tabulate_only

# Usage: Rscript cluster_by_taxon_p1.R [primer_pair] [infilename] [kingdom_column] [threads]
#
if [[ $aggregate_only == "F" ]]
then
	echo "Tabulating taxa..."
	echo "Running: Rscript cluster_by_taxon_p1.R $primer_pair $psfilename $kindgom_column $threads"
	#Rscript cluster_by_taxon_p1.R $primer_pair $psfilename $kingdom_column $threads
	Rscript cluster_by_taxon_p1_parallel.R $primer_pair $psfilename $kingdom_column $threads

fi


# Usage: Rscript cluster_by_taxon_p2.R [primer_pair] [infilename] [singleton_cutoff] [taxon_level] [rarefaction_level] [threads] [kingdom_column] [fasta_folder] [resume]
#
if [[ $tabulate_only == "F" ]]
then
	echo "Aggregating reads by taxon..."
	echo "Running: Rscript cluster_by_taxon_p2.R $primer_pair $psfilename $singleton_cutoff $taxonlevel $rarefaction_level $threads $kingdom_column $fasta_folder $resume"
	Rscript cluster_by_taxon_p2_parallel.R $primer_pair $psfilename $singleton_cutoff $taxonlevel $rarefaction_level $threads $kingdom_column $fasta_folder $resume
fi

# if you made it to here without an error, write to session log file
echo "aggregate_minimap_bytaxon.sh completed without error" >> $session_log_file
echo "" >> $session_log_file