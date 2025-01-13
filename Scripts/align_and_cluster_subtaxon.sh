#!/bin/bash
# February 6, 2024

	# Usage: sh align_and_cluster_subtaxon.sh [-t subtaxon_alignment_threads] [-m mafft_threads] [-p phylip_threads] [-c mothur_threads] [-P mothur_path] [-A alignment only] [-D distance only] [-S subtaxon clustering only] [-B post subtaxon bin seqs only] [-C consensus sequences only] [-W pairwise alignment only] [-R resume task] [-i input directory to look for fasta files] [-o output directory for mothur] [-L session log file] [-h display this help message]

	# Description: This script is a dynamic command that performs several tasks to ulimately cluster and generate consensus sequences:
		# i. generates mafft alignment-based OR mothur pairwise.dist-based distance matrices
				# i.a.1 mafft
				# i.a.2 phylip dna.dist
					# OR
				# i.b mothur pairwise.dist
		# ii. from distances in (i), "subtaxon" clustering with mothur cluster.classic [NOTE: need to implement directory specification]
		# iii. "bin seqs" with mothur bin.seqs at the chosen distance cutoff label
		# iv.  subsequent alignments of sequences within the subtaxon clusters and generation of consensus sequences
	
	# Options (note some options are currently applied for different specifications in multiple steps):
	
		# General Options:
			# -P mothur_path='~/mothur/mothur'
			# -i input_fasta_directory='taxon_cluster' 		- input directory for mothur to look for fasta files
			# -o output_mothur_directory='distance_matrices' 	- directory where mothur saves output files (may need to change depending on steps being run)
			# -A 												- only performs the first MAFFT alignment step (i.a.1)
			# -D												- only performs the phylip dnadist step (i.a.2)
			# -S												- only performs "subtaxon" clustering (ii) in mothur
			# -B												- only runs mothur bin.seqs
			# -C												- only does subtaxon alignments and consensus sequence calling
			# -W												- generates pairwise distance matrices in mothur with pairwise.dist
			# -L session_log_file=""							- specifies name of log file (default align_and_cluster_subtaxon.log)
	
		# MAFFT and phylip alignment and distance options (i.a)
			# -t subtaxon_alignment_threads=4	- number of instances of mafft to run
			# -m mafft_threads=4				- number of threads to run in mafft
			# -p phylip_threads=16			- number of instances of phylip to run
			# -R resume='F'					- pick up where left off (currently only implemented here)
		
		# mothur pairwise.dist Options
			# -c mothur_threads=16			- number of instance of mothur to run
		
		# Subtaxon clustering (mothur cluster.classic)
			# -t subtaxon_alignment_threads=4	- number of instances of mothur to run
			# [NOTE: need to implement directory specification]

	# Output files:
		# i.a.1	taxon_alignments/*.maf
		# i.a.2	distance_matrices/*.maf.dist
		# i.b		output_mothur_directory/*.phylip.dist
		# ii.		*.rabund [NOTE: need to implement directory specification]
				# *.sabund [NOTE: need to implement directory specification]
				# *.list [NOTE: need to implement directory specification]
		# iii.	*.an.[cutoff].fasta [NOTE: need to implement directory specification]
		# iv.		*.cons.denovo.python.fasta	- final seq file
				# *.cons.denovo.python.data	- data file
				# *.reformatted.fasta			- seq file
				# *.karuna.cons.ungap			- seq file
				# *.karuna.cons				- seq file
				# *.cons.data					- data file
# options
subtaxon_alignment_threads=4
mafft_threads=4
phylip_threads=16
mothur_threads=16
mothur_path='~/mothur/mothur'
input_fasta_directory='taxon_cluster'
output_mothur_directory='distance_matrices'
resume='F'
session_log_file=""

# internal 'skip' parameters
pairwise_alignment='F'
skip_alignment='F'
skip_distance='F'
skip_subtaxon_clustering='F'
skip_bin_seqs='F'
skip_consensus_sequences='F'

print_usage(){
	echo "Usage: sh align_and_cluster_subtaxon.sh [-t subtaxon_alignment_threads] [-m mafft_threads] [-p phylip_threads] [-c mothur_threads] [-P mothur_path] [-A alignment only] [-D distance only] [-S subtaxon clustering only] [-B post subtaxon bin seqs only] [-C consensus sequences only] [-W pairwise alignment only] [-R resume task] [-i input directory to look for fasta files] [-o output directory for mothur] [-L session log file] [-h display this help message] " >&2
}
# 1. Create ProgressBar function
# 1.1 Input is currentState($1) and totalState($2)
function ProgressBar {
# Process data
    let _progress=(${1}*100/${2}*100)/100
    let _done=(${_progress}*4)/10
    let _left=40-$_done
# Build progressbar string lengths
    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")

# 1.2 Build progressbar strings and print the ProgressBar line
# 1.2.1 Output example:                           
# 1.2.1.1 Progress : [########################################] 100%
	printf "Progress : [${_fill// /#}${_empty// /-}] ${_progress}%%"

}
while getopts 't:m:p:c:P:ADSBCWi:o:RL:h' OPTION; do
	case "$OPTION" in
		t)
			subtaxon_alignment_threads="$OPTARG"
			;;
		m) mafft_threads="$OPTARG" ;;
		p) phylip_threds="$OPTARG" ;;
		c) mothur_threads="$OPTARG" ;;
		P) mothur_path="$OPTARG" ;;
		A)
			skip_distance='T'
			skip_subtaxon_clustering='T'
			skip_bin_seqs='T'
			skip_consensus_sequences='T'
			;;
		D)
			skip_alignment='T'
			skip_subtaxon_clustering='T'
			skip_bin_seqs='T'
			skip_consensus_sequences='T'
			;;
		S)
			skip_alignment='T'
			skip_distance='T'
			skip_bin_seqs='T'
			skip_consensus_sequences='T'
			;;
		B)
			skip_alignment='T'
			skip_distance='T'
			skip_subtaxon_clustering='T'
			skip_consensus_sequences='T'
			;;
		C)
			skip_alignment='T'
			skip_distance='T'
			skip_subtaxon_clustering='T'
			skip_bin_seqs='T'
			;;
		W)
			pairwise_alignment='T'
			skip_alignment='T'
			skip_distance='T'
			skip_subtaxon_clustering='T'
			skip_bin_seqs='T'
			skip_consensus_sequences='T'
			;;
		i) input_fasta_directory="$OPTARG" ;;
		o) output_mothur_directory="$OPTARG" ;;
		R)
			resume='T'
			;;
		L) session_log_file="$OPTARG" ;;
		h)
			print_usage
			exit ;;
		*)
		exit 1 ;;
	esac
done

UPLINE=$(tput cuu1)
ERASELINE=$(tput el)

##### WRITE EVERYTHING TO AN OUTPUT FILE
# CURRENTLY DISABLED
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
file_name="align-and-cluster-subtaxon-log"
file_name_stamped=$file_name.$current_time.txt

#{

echo "number of worker threads (subtaxon_alignment_threads) = $subtaxon_alignment_threads"
echo "mafft_threads = $mafft_threads"

#{
#echo "number of worker threads (subtaxon_alignment_threads) = $subtaxon_alignment_threads"
#echo "mafft_threads = $mafft_threads"
#}&> "$file_name_stamped"

# if session log file name not supplied, create one (e.g., ITS54.log)
if [[ $session_log_file == "" ]]; then
	log_file='align_and_cluster_subtaxon.log'
fi

# if session log file by that name does not yet exist, create a file
if ! [[ -f $session_log_file ]]; then
	>$session_log_file
fi

# next output command and options called to the session log file:
echo $current_time: >> $session_log_file
echo "sh align_and_cluster_subtaxon.sh -t $subtaxon_alignment_threads -m $mafft_threads -p $phylip_threads -c $mothur_threads -P $mothur_path -R $resume task] -i $input_fasta_directory -o $output_mothur_directory -L $session_log_file" >> $session_log_file
echo "pairwise_alignment=$pairwise_alignment">> $session_log_file
echo "skip_alignment (mafft) =$skip_alignment">> $session_log_file
echo "skip_distance (mafft and phylip) =$skip_distance">> $session_log_file
echo "skip_subtaxon_clustering (mothur) =$skip_subtaxon_clustering">> $session_log_file
echo "skip_bin_seqs (mothur) =$skip_bin_seqs">> $session_log_file
echo "skip_consensus_sequences=$skip_consensus_sequences">> $session_log_file
echo "..." >> $session_log_file

cd cluster_by_taxon

current_path=$(pwd)

###################################################
# pairwise alignments in mothur for higher taxons #
###################################################

if [[ $pairwise_alignment == 'T' ]]
then
	printf "Running pairwise dist...\n\n"
	if [ ! -d $output_mothur_directory ]
	then
		mkdir $output_mothur_directory
	fi
	#cd $input_fasta_directory
	
	let COUNTER=0
	let _number=0
	_end=$(ls -1q cluster_*.fa | wc -l)
	
	for fasta_file in cluster_*.fa
	do
		echo "$UPLINE$ERASELINE"
		printf "mafft aligning $fasta_file ($_number of $_end)...\n"
		ProgressBar ${_number} ${_end}	
	
		clusternumber=${fasta_file%.fa}

		# RENAME IT .mothur.pairwise.dist
	
		$mothur_path "#set.dir(input=$current_path$input_fasta_directory, output=$current_path$output_mothur_directory);pairwise(fasta=$fasta_file,format=lt, cutoff=100, processors=$mothur_threads);cluster.classic(phylip=$clusternumber.phylip.dist, cutoff=100)" &
	
		((COUNTER=COUNTER+1))
		if (($COUNTER == $mothur_threads))
		then
			wait -n
			((COUNTER=COUNTER-1))
		fi
		((_number=_number+1))	
	done
fi

#############
# alignment #
#############

if [[ $skip_alignment == 'F' ]]
then

	printf "Making alignments...\n\n"
	#echo "Making alignments..." &> "$file_name_stamped"

	if [ ! -d "taxon_alignments" ]
	then
		mkdir taxon_alignments
	fi

	if [[ $resume == 'T' ]]
	then
		cd taxon_alignments
		empties=$(find . -size 0)
		if [[ ! $empties == '' ]]
		then
			rm $empties
		fi
		cd ..
	fi

	cd taxon_cluster

	let COUNTER=0
	let _number=0
	_end=$(ls -1q cluster_*.fa | wc -l)
	
	for fasta_file in cluster_*.fa
	do
		if [[ $resume == 'T' ]]
		then
			if [ -f "../taxon_alignments/$fasta_file.maf" ]
			then
				((_number=_number+1))
				continue
			fi
		fi
		#mafft --globalpair --thread 8 --op 0.5 --gop -0.5 --phylipout $fasta_file > ../taxon_alignments/$fasta_file.maf
		#echo "mafft aligning $fasta_file..." &> "$file_name_stamped"
		echo "$UPLINE$ERASELINE"
		printf "mafft aligning $fasta_file ($_number of $_end)...\n"
		ProgressBar ${_number} ${_end}
		#{
			ginsi --thread $mafft_threads --phylipout --quiet $fasta_file > ../taxon_alignments/$fasta_file.maf &
		#}&> "$file_name_stamped"
		((COUNTER=COUNTER+1))
		if (($COUNTER == $subtaxon_alignment_threads))
		then
			wait -n
			((COUNTER=COUNTER-1))
		fi
		((_number=_number+1))
	done

	cd ..

fi

#####################
# distance matrices #
#####################

if [[ $skip_distance == 'F' ]]
then

	echo "Running phylip..."

	if [ ! -d "distance_matrices" ]
	then
		mkdir distance_matrices
	fi

	if [[ $resume == 'T' ]]
	then
		cd distance_matrices
		empties=$(find . -size 0)
		if [[ ! $empties == '' ]]
		then
			rm $empties
		fi
		cd ..
	fi

	cd taxon_alignments

	let COUNTER=0
	let _number=0
	_end=$(ls -1q cluster_*.maf | wc -l)

	echo "blah" > outfile
	for maf_file in cluster_*.maf
	do
		if [[ $resume == 'T' ]]
		then
			if [ -f "../distance_matrices/$maf_file.dist" ]
			then
				((_number=_number+1))
				continue
			fi
		fi
		echo "$UPLINE$ERASELINE"
		printf "dnadist $maf_file...\n"
		ProgressBar ${_number} ${_end}
		
		log_phylip=$maf_file.phylip.log
		input_phylip=$maf_file.phylip.input
		
		echo "$maf_file" > $input_phylip
		echo "F" >> $input_phylip
		echo "../distance_matrices/$maf_file.dist" >> $input_phylip
		echo "D" >> $input_phylip
		echo "D" >> $input_phylip
		echo "L" >> $input_phylip
		echo "Y" >> $input_phylip
		../../dnadist < $input_phylip &> $log_phylip &
		((COUNTER=COUNTER+1))
		if (($COUNTER == $phylip_threads))
		then
			wait -n
			((COUNTER=COUNTER-1))
		fi
		((_number=_number+1))
	done

	cd ..

fi

##########################
# subtaxon clustering... #
##########################

if [[ $skip_subtaxon_clustering == 'F' ]]
then

	if [ ! -d $output_mothur_directory ]
	then
		mkdir $output_mothur_directory
	fi

	cd distance_matrices

	echo "Running mothur cluster.classic..."

	let COUNTER=0
	for dist_file in cluster_*.maf.dist
	do
		clusternumber=${dist_file%\.maf\.dist}
		
		# this produces a file called cluster_*.an.list
		# also produces files called .an.rabund and .an.sabund
		$mothur_path "#set.dir($current_path$input_fasta_directory);cluster.classic(phylip=./$dist_file)" &
		
		((COUNTER=COUNTER+1))
		if (($COUNTER == $subtaxon_alignment_threads))
		then
			wait -n
			((COUNTER=COUNTER-1))
		fi
	done

	cd ..

fi

################################################################
################################################################

################# 11/14/2023 ####################################

# here we will want to decide how we are going to split them up
# we are gonna want to analyze the clustering results... across the
# data to find the ideal cutoff...

####################################################################

# bin seqs

if [[ $skip_bin_seqs == 'F' ]]
then

	if [ ! -d "cluster_files" ]
	then
		mkdir cluster_files
	fi

	cd distance_matrices

	echo "Running mothur bin.seqs..."

	let COUNTER=0
	for dist_file in cluster_*.maf.dist
	do
		clusternumber=${dist_file%\.maf\.dist}

		# this produces a lot of files...
		# going to want to use the above to chose the optimal label parameter
		echo "set.dir(output=../cluster_files);bin.seqs(list=./$dist_file.an.list, fasta=../taxon_cluster/$clusternumber.fa)" | $mothur_path
		
		((COUNTER=COUNTER+1))
		if (($COUNTER == $phylip_threads))
		then
			wait -n
			((COUNTER=COUNTER-1))
		fi
	done

	cd ..

fi

###########################
# consensus sequences ... #
###########################



####### THE REST NEEDS TO BE EDITED ITS JUST BEEN PASTED FROM OLD SCRIPT
####### OR MOVED TO A NEW SCRIPT

# then get consensus sequences for each otu in each taxon
# loop through the taxa

# will need another script to chose or optimize clustering percent... 
# make an R plot of number of clusters vs. cutoff
if [[ $skip_consensus_sequences == 'F' ]]
then
	# need to break up OTUs
	Rscript ../metabarcoding_consensus_files.R # need to have a way to set the cutoff!
	otus=$(ls Otu*.fasta | wc | awk '{print $1}')
	# loop through the otus
	if (( $otus == 0 ))
	then
		cd ..
			rm -R ${ARRAY_FOLDERS[j]}
	else
		for otu_file in Otu*.fasta
		do
			# gonna have to move the file here and then move it back
			ginsi --clustalout --inputorder --thread 64 --threadtb 20 --threadit 20 --quiet ${otu_file%.*}.concat.fasta > all_filt_denovo.maf
			
			# need to edit this so it outputs files with different names
			python cons.denovo.py #Otu*.maf 
			
			mv all_filt_denovo.maf $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.maf
			
			# also edit this to be compatible with the above
			Rscript consensus_seq_and_score.denovo.R $otu_file ${otu_file%.*}	 
			
			# change all the folder name stuff
			mv all_filt_denovo.maf.karuna.cons $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.maf.karuna.cons
			mv all_filt_denovo.maf.cons.data $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.maf.cons.data
			mv all_filt_denovo.maf.karuna.cons.ungap $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.maf.karuna.cons.ungap
			mv ${ARRAY_FOLDERS[j]}.cons.denovo.python.fasta $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.cons.denovo.python.fasta
			mv ${ARRAY_FOLDERS[j]}.cons.denovo.python.data $folder.${ARRAY_FOLDERS[j]}.${otu_file%.*}.cons.denovo.python.data
		done
	fi

fi

#} &> "$file_name_stamped"


# if you made it to here without an error, write to session log file
echo "align_and_cluster_subtaxon.sh completed without error" >> $session_log_file
echo "" >> $session_log_file


# then we will have to
#
# 1. Use TBAS (API?) to ID consensus sequences from LSU-ITS
#
##### advantage of this; TBAS trees can be added to final phyloseq object to generate unifrac distances!
#
#    XXXX BLAST- assign identities to the sequences (could be tricky) with NCBI
#
# 1a. for full ITS/simple ITS without LSU, could still use UNITE database
#
# 2. Create a new database with the consensus sequences (and putative identites)
#
# 3. is this necessary - ??Re- alignment search sequences (that are identified to the level at which we are limiting our clusting) using MiniMap2
#    ** COULD ALSO just use the seqid numbers to get the abundances of the different ones in each barcode
#
# 4. Compare pre- to post-assembly distances (what is the relationship?)

# 5. Try to add at Family/Order/higher levels!!