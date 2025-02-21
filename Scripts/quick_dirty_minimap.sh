# in this one we just do minimap alignments instead of the clustering approach
# tell the thing where the database is

#	Usage: sh quick_dirty_minimap.sh [-p primer_pair] [-q quality_cutoff] [-d database] [-i infile] [-b barcode_threads] [-t minimap_threads] [-m minlen] [-s skip_minimap] [-S skip_samtools]  [-L session_log_file]

#	Description: This script does the following
#		i.		runs minimap2 (a mapping alingment search) for each barcode to the database of choice and produces samtools output
#		ii.		creates the samtools output summary for each barcode
#		iii.	runs process_minimap2.R to convert that into a OTU table
#				the R script parses the samtools output to make it readable in R,
#				particularly in make_phyloseq.R

#	Options (* required) (= default)

#		General Options
#			-p (=primer_pair) 						- name of primer pair as it appears in primer_seqs.sh (Default 'ITS54')
#			-q (=quality_cutoff)					- the quality cutoff (as phred score) to filter initial sequences ((Default 10, corresponds to 90%)
#			-m minlen=300							- minimum length of seqs to filter initial sequences
#			-s skip_minimap - (Default 'F') can be used to skip the minimap step in combination with -S to skip the samtools summary of minimap step - note that skipping minimap but not samtools makes little sense unless you wanted to just do one step at a time; but the summary section will still run at the end to gather minimap output.
#			-S skip_samtools - (Default 'F') can be used in combination with -s to skip minimap and processing of minimap results. effectively, this function only gathers and summarizes minimap results.
#			-L session_log_file="[-p primer_pair].log" 
#				name of log file where you are keeping track of all the commands
#				you are running and with which parameters (automatically supplies
#				one based on primer pair)
	
#		MiniMap2 Options
#			-d* (=database)							- the name and location of reference database file for taxonomic assignments (formatted as a fasta file)
#			-i infile='all_filt_reorient.fasta'		- name of input fasta file in each barcode folder
#			-b barcode_threads=4					- number of minimaps to run independently (total cpus = b * t)
#			-t minimap_threads=4					- number of minimap threads to tell minimap run in each minimap instance
#			-P minimap_path="/vol3/home/ec2-user/minimap2-2.24_x64-linux"

#	Output files:
#		barcodeXX/primer_pair/all_filt_minimap.sam (minimap2 output)
#		barcodeXX/primer_pair/all_filt.samview.tsv (samtools view output)
#		barcodeXX/primer_pair/relabund_phred_q10.csv (database accession - like an OTU - raw abundance table)

# first set default values for parameters
primer_pair='ITS54'
quality_cutoff=10         #90%
database='' #'/vol3/home/ec2-user/ref_dbs/UNITE.fasta'
infile='all_filt_reorient.fasta'
barcode_threads=4
minimap_threads=4
minlen=300
minimap_path='minimap2' #"/vol3/home/ec2-user/minimap2-2.24_x64-linux"
skip_minimap='F'
skip_samtools='F'
session_log_file=""

# define a function to call if user asks linux to tell how to use the script
print_usage(){
	echo "Usage: sh quick_dirty_minimap.sh [-p primer_pair] [-q quality_cutoff] [-d database] [-i infile] [-b barcode_threads] [-t minimap_threads] [-m minlen] [-P minimap_path] [-s skip_minimap] [-S skip_samtools] [-L session_log_file]" >&2
}

# get the options from the tags
while getopts 'p:q:d:i:b:t:m:P:s:S:L:h' OPTION; do
	case "$OPTION" in
		p)
			primer_pair="$OPTARG"
			;;
		q) quality_cutoff="$OPTARG" ;;
		d) database="$OPTARG" ;;
		i) infile="$OPTARG" ;;
		b) barcode_threads="$OPTARG" ;;
		t) minimap_threads="$OPTARG" ;;
		m) minlen="$OPTARG" ;;
		P) minimap_path="$OPTARG" ;;
		s) skip_minimap="$OPTARG" ;;
		S) skip_samtools="$OPTARG" ;;
		L) session_log_file="$OPTARG" ;;
		h)
			print_usage
			exit ;;
		*)
		exit 1 ;;
	esac
done

# create a log file JUST FOR THIS CALL OF TRIM_AND_SORT with a unique time stamp
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
file_name="quick-dirty-minimap-log"
file_name_stamped=$file_name.$current_time.log

# if session log file name not supplied, create one (e.g., ITS54.log)
if [[ $session_log_file == "" ]]; then
	session_log_file=$file_name_stamped
fi

# if session log file by that name does not yet exist, create a file
if ! [[ -f $session_log_file ]]; then
	>$session_log_file
fi

# from here on, all output is rerouted to the log file unless otherwise specified
{
# first output command call to the log file just for this call
echo "running sh quick_dirty_minimap.sh -p $primer_pair -q $quality_cutoff -d $database -i $infile -b $barcode_threads -t $minimap_threads -m $minlen -P $minimap_path -s $skip_minimap -S $skip_samtools -L $session_log_file"

((total_threads = minimap_threads * barcode_threads))



# next output command and options to the session log file:
echo $current_time: >> $session_log_file
echo "sh quick_dirty_minimap.sh -p $primer_pair -q $quality_cutoff -d $database -i $infile -b $barcode_threads -t $minimap_threads -m $minlen -P $minimap_path -s $skip_minimap -S $skip_samtools -L $session_log_file" >> $session_log_file
echo "..." >> $session_log_file

# check to see if we are skipping minimap
# runs minimap2 (a mapping alingment search) for each barcode to the database of choice and produces samtools output
if [[ $skip_minimap == "F" ]]
then
	if [ -f $database ]; then
		let COUNTER=0
		
		# cycle through folders
		for folder in ./barcode*/
		do
			if [ -d "$folder/$primer_pair" ]
			then
				echo "minimap processing folder $folder"
				cd $folder
				cd $primer_pair
	#			cp ../../$database .

				# runs minimap2 (a mapping alingment search) for each barcode to the database of choice and produces samtools output
				$minimap_path/minimap2 -ax map-ont -t $minimap_threads --secondary=no $database $infile > all_filt_minimap.sam &
				((COUNTER=COUNTER+1))
				if (($COUNTER == $barcode_threads))
				then
					wait
					let COUNTER=0
				fi
		#		rm $database
				cd ../..
			fi
		done
		wait
	else
		echo "No reference database"
	fi
fi

# check to see if we are skipping samtools
# creates the samtools output summary for each barcode
if [[ $skip_samtools == "F" ]]
then
	let COUNTER=0
	
	# cycle through folders
	for folder in ./barcode*/
	do
		if [ -d "$folder/$primer_pair" ]
		then
			echo "samtools processing folder $folder"
			cd $folder
			cd $primer_pair
			
			# creates the samtools output summary for each barcode
			samtools view all_filt_minimap.sam -o all_filt.samview.tsv &
			((COUNTER=COUNTER+1))
			if (($COUNTER == $total_threads))
			then
				wait
				let COUNTER=0
			fi
			cd ../..
		fi
	done
	wait
fi

# check to see if we are skipping samtools
# runs process_minimap2.R to convert that into a OTU table
let COUNTER=0
for folder in ./barcode*/
do
	if [ -d "$folder/$primer_pair" ]
	then
		echo "summarizing minimap output in folder $folder"
		cd $folder
		cd $primer_pair
		
		# runs process_minimap2.R to convert that into a OTU table
		Rscript ../../process_minimap2.R $quality_cutoff $minlen &
		((COUNTER=COUNTER+1))
		if (($COUNTER == $total_threads))
		then
			wait
			let COUNTER=0
		fi
		cd ../..
	fi
done

} &> "$session_log_file"

# if you made it to here without an error, write to session log file
echo "quick_dirty_minimap.sh completed without error" >> $session_log_file
echo "" >> $session_log_file