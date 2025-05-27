# sh make_joined_fasta.sh

#!/bin/bash
#	Usage: sh make_joined_fasta.sh [-p primer_pair] [-t taxidfilename] [-s taxideseqfilename] [-a aggregated_fasta_file] [-L session_log_file]

#default parameter values
primer_pair='ITS54'
taxidfilename='cluster_by_taxon/taxids.tsv'
taxideseqfilename='cluster_by_taxon/taxid_seqid.tsv'
aggregated_fasta_file='all_sequences.fasta'

print_usage(){
	echo "Usage: sh aggregate_minimap_bytaxon.sh [-p primer_pair] [-t taxidfilename] [-s taxidseqfilename] [-a aggregated_fasta_file] [-L session_log_file]"
}
while getopts 'p:t:s:a:L:h' OPTION; do
	case "$OPTION" in
		p)
			primer_pair="$OPTARG"
			;;
		t) taxidfilename="$OPTARG" ;;
		s) taxidseqfilename="$OPTARG" ;;
		a) aggregated_fasta_file="$OPTARG" ;;
		L) session_log_file="$OPTARG" ;;
		h)
			print_usage
			exit ;;
		*)
		exit 1 ;;
	esac
done

current_time=$(date "+%Y.%m.%d-%H.%M.%S")
file_name="make-joined-fasta-log"
file_name_stamped=$file_name.$current_time.txt

# if session log file name not supplied, create one (e.g., ITS54.log)
if [[ $session_log_file == "" ]]; then
	log_file=$primer_pair.log
fi

# if session log file by that name does not yet exist, create a file
if ! [[ -f $session_log_file ]]; then
	>$session_log_file
fi

echo "running make_joined_fasta.sh -p $primer_pair -t $taxidfilename -s $taxidseqfilename -a $aggregated_fasta_file -L $session_log_file"
printf "running make_joined_fasta.sh -p $primer_pair -t $taxidfilename -s $taxidseqfilename -a $aggregated_fasta_file -L $session_log_file"

# next output command and options to the session log file:
echo $current_time: >> $session_log_file
echo "sh make_joined_fasta.sh -p $primer_pair -t $taxidfilename -s $taxidseqfilename -a $aggregated_fasta_file -L $session_log_file"
echo "..." >> $session_log_file

echo "Using $primer_pair"
echo "Using $taxidfilename and $taxidseqfilename"

cat barcode*/ITS54/all_filt_reorient.fasta > $aggregated_fasta_file

Rscript make_joined_fasta.R $primer_pair $taxidfilename $taxidseqfilename $aggregated_fasta_file



# if you made it to here without an error, write to session log file
echo "aggregate_minimap_bytaxon.sh completed without error" >> $session_log_file
echo "" >> $session_log_file