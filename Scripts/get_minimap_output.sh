#This reads in the minimap csv output from each folder and copies them all to a named file in the folder minimap_primer_pair

#	Usage: sh get_minimap_output.sh [-p primer_pair] [-L session_log_file]
	
#	Description: This reads in the minimap csv output from each folder and copies them all to a named file in the folder minimap_primer_pair
	
#	Options (* required)
#		-p* (=primer_pair) - name of primer pair as it appears in primer_seqs.sh
#		-L session_log_file="[-p primer_pair].log" 
#			name of log file where you are keeping track of all the commands
#			you are running and with which parameters (automatically supplies
#			one based on primer pair)
		
#	Output: ./minimap_primer_pair/*.barcodeXX.csv

primer_pair=''
session_log_file=''

print_usage(){
	echo "Usage: sh get_minimap_output.sh [-p primer_pair] [-L session_log_file]" >&2
}

# print options
while getopts 'p:L:h' OPTION; do
	case "$OPTION" in
		p) primer_pair="$OPTARG" ;;
		L) session_log_file="$OPTARG" ;;
		h)
			print_usage
			exit ;;
		*)
		exit 1 ;;
	esac
done

# check for primer name
if (( $primer_pair == '' )); then 
	echo "You must specify -p primer_pair folder name..."
	exit 1
fi

# if session log file name not supplied, create one (e.g., ITS54.log)
if (( $session_log_file == '' )); then
	session_log_file=get_minimap_output.$primer_pair.log
fi

# if session log file by that name does not yet exist, create a file
if ! [[ -f $session_log_file ]]; then
	>$session_log_file
fi

# next output command and options to the session log file:
echo $current_time: >> $session_log_file
echo "sh get_minimap_output.sh -p $primer_pair -L $session_log_file">> $session_log_file
echo "..." >> $session_log_file

# make the directory (will throw error if directory already exists)
mkdir minimap_$primer_pair

# cycle through folders and pull and rename the csv files
for folder in barcode*/
do
	cd $folder
#	echo $folder
	ff=$(echo $folder | tr -dc '0-9')
#	echo $ff
	if [ -d $primer_pair ]
	then
		cd $primer_pair
		for file in *.csv
		do
			if [ -f "$file" ]
			then
				fff=${file%.*}
				cp $file ../../minimap_$primer_pair/$fff.barcode$ff.csv
				#echo "$fff.barcode$ff.csv"
			fi
		done
		cd ..
	fi
	cd ..
done

#cp barcode*/$primer_pair/*.csv minimap_$primer_pair

# if you made it to here without an error, write to session log file
echo "get_minimap_output.sh completed without error" >> $session_log_file
echo "NOTE: must run Rscript make_phyloseq.R before proceeding to aggregate_minimap_bytaxon.sh">> $session_log_file
echo "" >> $session_log_file