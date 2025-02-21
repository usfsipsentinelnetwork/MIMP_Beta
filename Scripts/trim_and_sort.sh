#!/bin/bash
#
## This script is intended for the processing of amplicons from
## MinIon runs - first step: TRIM AND SORT
## needed in directory:
## ReverseComplement.pl
## needed loaded/installed into space/module:
## NanoFilt
## NanoPlot
## cutadapt

#	Usage: sh trim_and_sort.sh [-p primer_pair] [-q quality_cutoff] [-m min_length] [-M max_length] [-c head_crop] [-a adapter_error] [-o adapter_overlap] [-N nanoplot_all(T/F)] [-n nanoplot_each(T/F)] [-x cpu_cores] [-s skip_to_cut_adapt] [-L session_log_file]

#	Description: This script does the following
#		i.	runs nanoplot (summary and visualization of quality and length, distribution etc. after basecalling)
#		ii.	runs nanofilt (filters by q score, sequence length, and trims some reads like adapters from the reads)
#		iii.	converts file to fasta format and shortens names of sequences in fasta sample headers for downstream compatibility
#		iv.	makes a new folder for the results for the primer set being used
#		v. 	runs cutadapt to trim forward and reverse primer seqs and reorients them and contatenates all seqs from each barcode

#	Options (and default values)

#		General
#			-p primer_pair='ITS1F4'				- the primer pair you want to use
#			-x skip_to_cut_adapt='F'			- whether to skip nanoplot & nanofilt and go straight to cutadapt
#			-L session_log_file="[-p primer_pair].log" 
#				name of log file where you are keeping track of all the commands
#				you are running and with which parameters (automatically supplies
#				one based on primer pair)

#		NanoFilt Options
#			-q quality_cutoff=15	- the quality cutoff for nanofilt
#			-m min_length=200 	- the minimum length for nanofilt
#			-M max_length=2000	- the maximum read length for nanofilt
#			-c head_crop=50		- the number of leading bases for nanofilt to trim

#		CutAdapt Options
#			-e adapter_error=.1   	- adapter error rate
#			-o adapter_overlap=15	- the amount of overlap allowed
#			-x cpu_cores=8			- number of cores to use for cut adapt

#		NanoPlot Options
#			-N nanoplot_all='T'	- whether to make a nanoplot for all seqs combined
#			-n nanoplot_each='T'	- whether to make a nanoplot for each barcode

#	Output files:
#		barcodeXX/nanoplot files
#		barcodeXX/all_filt_concatenated.fasta (before cutadapt)
#		barcodeXX/primer_pair/all_filt_reorient.fasta

# first set default values for parameters
primer_pair='ITS1F4'
quality_cutoff=15
min_length=200
max_length=2000
head_crop=50
adapter_error=.1
adapter_overlap=15 # double check this
nanoplot_all='T'
nanoplot_each='T'
cpu_cores=8
skip_to_cut_adapt='F'
session_log_file=""

# define a function to call if user asks linux to tell how to use the script
print_usage(){
	echo "Usage: sh trim_and_sort.sh [-p primer_pair] [-q quality_cutoff] [-m min_length] [-M max_length] [-c head_crop] [-a adapter_error] [-o adapter_overlap] [-N nanoplot_all(T/F)] [-n nanoplot_each(T/F)] [-x cpu_cores] [-s skip_to_cut_adapt] [-L session_log_file]" >&2
}

# get the options from the tags
while getopts 'p:q:m:M:c:e:o:N:n:x:s:L:h' OPTION; do
	case "$OPTION" in
		p)
			primer_pair="$OPTARG"
			;;
		q) quality_cutoff="$OPTARG" ;;
		m) min_length="$OPTARG" ;;
		M) max_length="$OPTARG" ;;
		c) head_crop="$OPTARG" ;;
		e) adapter_error="$OPTARG" ;;
		o) adapter_overlap="$OPTARG" ;;
		N) nanoplot_all="$OPTARG" ;;
		n) nanoplot_each="$OPTARG" ;;
		x) cpu_cores="$OPTARG" ;;
		s) skip_to_cut_adapt="$OPTARG" ;;
		L) session_log_file="$OPTARG" ;;
		h)
			print_usage
			exit ;;
		*)
		exit 1 ;;
	esac
done

# read in primer file
. ./primer_seqs.sh

# create an array that has the folder names in it
index=$(echo ${ARRAY_FOLDERS[@]/$primer_pair//} | cut -d/ -f1 | wc -w | tr -d ' ')

# create a log file JUST FOR THIS CALL OF TRIM_AND_SORT with a unique time stamp
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
file_name="trim_and_sort-log"
file_name_stamped=$file_name.$current_time.log

# from here on, all output is rerouted to the log file unless otherwise specified
{
# first output command call to the log file just for this call
echo "running trim_and_sort.sh -p $primer_pair -q $quality_cutoff -m $min_length -M $max_length -c $head_crop -e $adapter_error -o $adapter_overlap -N $nanoplot_all -n $nanoplot_each -x $cpu_cores -s $skip_to_cut_adapt -L $session_log_file"

# first need to check to make sure the primer name used is in the file
if [[ $index == -1 ]]; then
	echo "no primerset named $primer_pair"
	exit 1
fi

# if session log file name not supplied, create one (e.g., ITS54.log)
if [[ $session_log_file == "" ]]; then
	session_log_file=$file_name_stamped
fi

# if session log file by that name does not yet exist, create a file
if ! [[ -f $session_log_file ]]; then
	>$session_log_file
fi

# next output command and options to the session log file:
echo $current_time: >> $session_log_file
echo "trim_and_sort.sh -p $primer_pair -q $quality_cutoff -m $min_length -M $max_length -c $head_crop -e $adapter_error -o $adapter_overlap -N $nanoplot_all -n $nanoplot_each -x $cpu_cores -s $skip_to_cut_adapt -L $session_log_file" >> $session_log_file
echo "..." >> $session_log_file

# runs nanoplot FOR ALL THE READS (summary and visualization of quality and length, distribution etc. after basecalling)
if [[ $nanoplot_all == 'T' ]]; then
		# b) generate quality plots
		if [ -d nanoplot ]; then
			rm -r nanoplot
		fi
		mkdir nanoplot
		NanoPlot --fastq barcode*/*.fastq --plots dot --outdir nanoplot
fi

# cycle through the barcode folders
#for folder in "${barcodes[@]}"
for folder in ./barcode*/
do

	# check to see if folder exists first
	if [ -d "$folder" ]
	then

		cd $folder
		
		if [[ $skip_to_cut_adapt == "F" ]]; then
			
			count=`ls all*fas* 2>/dev/null | wc -l`
			if (( $count > 0 )); then
				rm all*fas*
			fi 
			
			count=`ls *.fastq 2>/dev/null | wc -l`
			if (( $count > 0 )); then
				cat *.fastq > all.fastq
				
				# runs nanoplot FOR EACH BARCODE (summary and visualization of quality and length, distribution etc. after basecalling)
				if [[ $nanoplot_each == 'T' ]]; then
				
					# b) generate quality plots
					echo "Running nanoplot on $folder for ${ARRAY_FOLDERS[index]}"
					if [ -d nanoplot ]; then
						rm -r nanoplot
					fi
					mkdir nanoplot
					NanoPlot --fastq all.fastq --plots dot --outdir nanoplot
				fi
				
				# runs nanofilt (filters by q score, sequence length, and trims some reads like adapters from the reads)
				# Q, minimum length bp, trim leading bases (adapter)
				echo "Running nanofilt on $folder for ${ARRAY_FOLDERS[index]}"
				NanoFilt -q $quality_cutoff -l $min_length --maxlength $max_length --headcrop $head_crop all.fastq > all_filt.fastq
				
				#d) convert to .fasta format
				echo "Converting $folder to fasta format for ${ARRAY_FOLDERS[index]}"
				sed -n '1~4s/^@/>/p;2~4p' all_filt.fastq > all_filt.fasta
				rm all.fastq all_filt.fastq
				
				#e) shorten names of sequences to just 8 characters
				sed 's/.//10g; n' all_filt.fasta > all_filt_concatenated.fasta
			else
				echo "No fastq files in $folder. make sure they are decompressed."
			fi
		fi
		echo "Running cutadapt on $folder for ${ARRAY_FOLDERS[index]}"
		
		if [ -f "all.fasta" ]; then

			DIRECTORY=${ARRAY_FOLDERS[index]}
			
			if [ -d $DIRECTORY ]; then
				echo "Folder already exists for primer set $DIRECTORY in $folder. removing..."
				rm -Rf ${ARRAY_FOLDERS[index]}
				mkdir ${ARRAY_FOLDERS[index]}
			else
				echo "Making folder for $DIRECTORY in $folder..."
				mkdir ${ARRAY_FOLDERS[index]}
			fi

			cd ${ARRAY_FOLDERS[index]}
			
			# runs cutadapt to trim forward and reverse primer seqs and reorients them
			cutadapt --rc -j $cpu_cores -g ${PRIMERS_F[index]} -o forward_seqs.fasta -m ${PRODUCT_LENGTH_MIN[index]} -M ${PRODUCT_LENGTH_MAX[index]} --overlap $adapter_overlap -e $adapter_error ../all_filt_concatenated.fasta --untrimmed-output=all_filt_noF.fasta
			cutadapt --rc -j $cpu_cores -a ${PRIMERS_RRC[index]} -o reverse_seqs.fasta -m ${PRODUCT_LENGTH_MIN[index]} -M ${PRODUCT_LENGTH_MAX[index]} --overlap $adapter_overlap -e $adapter_error all_filt_noF.fasta --discard-untrimmed
			
			# contatenates all seqs from each barcode
			cat forward_seqs.fasta reverse_seqs.fasta > all_filt_reorient.fasta
			myfilesize1=$(wc -l "all_filt_reorient.fasta" | awk '{print $1}')
			cd ..
			
			# remove file with less than 3 reads for the whole barcode
			if (( $myfilesize1 < 6 ))
			then
				rm -Rf ${ARRAY_FOLDERS[index]}
			fi
		fi
		cd ..
	fi
done
} > "$file_name_stamped"

# if you made it to here without an error, write to session log file
echo "trim_and_sort.sh completed without error" >> $session_log_file
echo "" >> $session_log_file