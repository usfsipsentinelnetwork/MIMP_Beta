LAST UPDATED FEB 6 2024

++++++++++++++++++++++++++++++++++++++++++++
++                                        ++
++   MIMP: MinIon Metabarcoding Pipeline  ++
++   ===================================  ++
++                                        ++
++      (C) Geoffrey Williams 2024        ++
++                                        ++
++++++++++++++++++++++++++++++++++++++++++++

+------------+
|  I. About  |
+------------+

MIMP (MinIon Metabarcoding Pipeline) features executable shell scripts to which one can pass arguments to augment a basic metabarcoding pipeline for processing multiplexed amplicon data from ONT MinIon, after basecalling. Currently the implementation is simple. MIMP supports two metabarcode amplicon pipelines a) "Quick & Dirty" and b) "De Novo" which respectively rely a) solely on alignment to the UNITE database or b) use preliminary UNITE alignments to inform within-taxon clustering (based on rarefied subsample of sequences) and re-alignment to resulting consensus sequences which can be BLASTed to a custom database (ie NCBI), for long-read amplicon metabarcoding. MIMP also supports c) processing amplicon data from single-sample libraries (MIMP "Sanger"). 

+--------------------+
|  II. Dependencies  |
+--------------------+

cutadapt
minimap2
samtools
python
R (and various R packages)
	dplyr (NOTE: early versions needed dplyr-devel but should not any more)
	tidyr
	phyloseq
	progress
	Biostrings
	foreach
	doSNOW
mafft
mothur

+------------------+
|  III. Workflows  |
+------------------+

a) MIMP "Quick & Dirty"
-----------------------

In order to run this pipeline, you should have your output from basecalling arranged into folders "barcode01" or "barcode01" through "barcode**" nad the following scripts in the folders, 

MIMP "Quick & Dirty": output OTU table/ phyloseq object solely on alignment to the UNITE database.

0. primer_seqs.sh

	Usage: NA

	Description: This is not a command to run; this file contains all the information about primers and primer sequences that you may want to select from using the command line call; this script can be edited to add additional information in the same format that it is set up for GNU BASH langauge.

	Each primer set (F and R) is given a number; here is the format used to specify seq, rev comp seq, and length min and max parameters:

	# ITS54
	# ITS5 GGAAGTAAAAGTCGTAACAAGG
	# ITS5(revcomp) CCTTGTTACGACTTTTACTTCC
	ARRAY_FOLDERS[4]="ITS54"
	PRIMERS_F[4]="GGAAGTAAAAGTCGTAACAAGG"
	PRIMERS_FRC[4]="CCTTGTTACGACTTTTACTTCC"
	PRIMERS_R[4]="TCCTCCGCTTATTGATATGC"
	PRIMERS_RRC[4]="GCATATCAATAAGCGGAGGA"
	PRODUCT_LENGTH_MIN[4]="300"
	PRODUCT_LENGTH_MAX[4]="1200"

1. trim_and_sort.sh

	Usage: sh trim_and_sort.sh [-p primer_pair] [-q quality_cutoff] [-m min_length] [-M max_length] [-c head_crop] [-a adapter_error] [-o adapter_overlap] [-N nanoplot_all(T/F)] [-n nanoplot_each(T/F)] [-x cpu_cores] [-s skip_to_cut_adapt] [-L session_log_file]

	Description: This script does the following
		i.	runs nanoplot (summary and visualization of quality and length, distribution etc. after basecalling)
		ii.	runs nanofilt (filters by q score, sequence length, and trims some reads like adapters from the reads)
		iii.	converts file to fasta format and shortens names of sequences in fasta sample headers for downstream compatibility
		iv.	makes a new folder for the results for the primer set being used
		v. 	runs cutadapt to trim forward and reverse primer seqs and reorients them and contatenates all seqs from each barcode

	Options (and default values)

		General
			-p primer_pair='ITS1F4'				- the primer pair you want to use
			-x skip_to_cut_adapt='F'			- whether to skip nanoplot & nanofilt and go straight to cutadapt
			-L session_log_file="[-p primer_pair].log" 
				name of log file where you are keeping track of all the commands
				you are running and with which parameters (automatically supplies
				one based on primer pair)

		NanoFilt Options
			-q quality_cutoff=15	- the quality cutoff for nanofilt
			-m min_length=200 	- the minimum length for nanofilt
			-M max_length=2000	- the maximum read length for nanofilt
			-c head_crop=50		- the number of leading bases for nanofilt to trim

		CutAdapt Options
			-e adapter_error=.1   	- adapter error rate
			-o adapter_overlap=15	- the amount of overlap allowed
			-x cpu_cores=8			- number of cores to use for cut adapt

		NanoPlot Options
			-N nanoplot_all='T'	- whether to make a nanoplot for all seqs combined
			-n nanoplot_each='T'	- whether to make a nanoplot for each barcode

	Output files:
		barcodeXX/nanoplot files
		barcodeXX/all_filt_concatenated.fasta (before cutadapt)
		barcodeXX/primer_pair/all_filt_reorient.fasta

2. quick_dirty_minimap.sh

	Usage: sh quick_dirty_minimap.sh [-p primer_pair] [-q quality_cutoff] [-d database] [-i infile] [-b barcode_threads] [-t minimap_threads] [-m minlen] [-s skip_minimap] [-S skip_samtools]  [-L session_log_file]

	Description: This script does the following
		i.		runs minimap2 (a mapping alingment search) for each barcode to the database of choice and produces samtools output
		ii.		creates the samtools output summary for each barcode
		iii.	runs process_minimap2.R to convert that into a OTU table
				the R script parses the samtools output to make it readable in R,
				particularly in make_phyloseq.R

	Options (* required) (= default)

		General Options
			-p (=primer_pair) 						- name of primer pair as it appears in primer_seqs.sh (Default 'ITS54')
			-q (=quality_cutoff)					- the quality cutoff (as phred score) to filter initial sequences ((Default 10, corresponds to 90%)
			-m minlen=300							- minimum length of seqs to filter initial sequences
			-s skip_minimap - (Default 'F') can be used to skip the minimap step in combination with -S to skip the samtools summary of minimap step - note that skipping minimap but not samtools makes little sense unless you wanted to just do one step at a time; but the summary section will still run at the end to gather minimap output.
			-S skip_samtools - (Default 'F') can be used in combination with -s to skip minimap and processing of minimap results. effectively, this function only gathers and summarizes minimap results.
			-L session_log_file="[-p primer_pair].log" 
				name of log file where you are keeping track of all the commands
				you are running and with which parameters (automatically supplies
				one based on primer pair)
	
		MiniMap2 Options
			-d* (=database)							- the name and location of reference database file for taxonomic assignments (formatted as a fasta file)
			-i infile='all_filt_reorient.fasta'		- name of input fasta file in each barcode folder
			-b barcode_threads=4					- number of minimaps to run independently (total cpus = b * t)
			-t minimap_threads=4					- number of minimap threads to tell minimap run in each minimap instance
			-P minimap_path="/vol3/home/ec2-user/minimap2-2.24_x64-linux"

	Output files:
		barcodeXX/primer_pair/all_filt_minimap.sam (minimap2 output)
		barcodeXX/primer_pair/all_filt.samview.tsv (samtools view output)
		barcodeXX/primer_pair/relabund_phred_q10.csv (database accession - like an OTU - raw abundance table)

3. get_minimap_output.sh

	Usage: sh get_minimap_output.sh [-p primer_pair] [-L session_log_file]
	
	Description: This reads in the minimap csv output from each folder and copies them all to a named file in the folder minimap_primer_pair
	
	Options (* required)
		-p* (=primer_pair) - name of primer pair as it appears in primer_seqs.sh
		-L session_log_file="[-p primer_pair].log" 
			name of log file where you are keeping track of all the commands
			you are running and with which parameters (automatically supplies
			one based on primer pair)
		
	Output: ./minimap_primer_pair/*.barcodeXX.csv

4. make_phyloseq.R ** required for future steps, but does not write to session log file

	Usage: Rscript make_phyloseq.R [primer_pair*] [db*] [outfilename] [silva_path]

	Description: Uses preliminary quick and dirty minimap2 to create an OTU and taxon table to make a phyloseq object

	Options (* required)
		1 [primer_pair] - the primer pair name
		2 [db]			- the database (currently supports =UNITE and =SILVA)
		3 [outfilename] - name of .RData output file with phyloseq object
		4 [silva_path]  - ***required if db=SILVA

	Output: outfilename.RData

b) MIMP "De Novo"
-----------------

MIMP "De Novo": use preliminary UNITE alignments to inform within-taxon clustering (based on rarefied subsample of sequences) and re-alignment to resulting consensus sequences which can be BLASTed to a custom database (ie NCBI), for long-read amplicon metabarcoding

5. aggregate_minimap_bytaxon.sh

	Usage: sh aggregate_minimap_bytaxon.sh [-p primer_pair] [-f phyloseq filename] [-l taxon level for clustering] [-h print this help message] [-s singleton_cutoff] [-r rarefaction_level] [-k kingdom_column] [-t threads] [-T tabulate only] [-A aggregate only] [-R resume] [-F fasta_folder]

	Description: Uses output from preliminary quick and dirty minimap2 to do run the following R scripts
		(supplying the options given to it to the R scripts); descriptions and specs for R scripts below

		5a. cluster_by_taxon_p1_parallel.R
			Usage: Rscript cluster_by_taxon_p1.R [primer_pair] [infilename] [kingdom_column] [threads]	
			Description: Tabulates seqids by taxon and outputs them to a file
			Options:
				1 [primer_pair] - the primer pair name (current default 'ITS1FLR3')
				2 [infilename]	- .RData file containing phyloseq object from quick_and_dirty and make_phyloseq
										Default: paste("minimap_",primer_pair,"/Phyloseq_Outfile_MIMP_",primer_pair,".RData",sep="")
				3 [kingdom_column] 	- informs the script which column of the in file tax table contains kingdom information, for compatibility with samtools/minimap2 output from various databases
				4 [threads]  		- number of threads to run	
			Output: cluster_by_taxon/taxid_seqid.tsv
		
		5b. cluster_by_taxon_p2_parallel.R	
			Usage: Rscript cluster_by_taxon_p2.R [primer_pair] [infilename] [singleton_cutoff] [taxon_level] [rarefaction_level] [threads] [kingdom_column] [fasta_folder] [resume]	
			Description: Uses output from p2 to aggregate sequences at the desired taxonomic and rarefaction levels
			Options:
				1 [primer_pair] 		- the primer pair name (current default 'ITS1FLR3')
				2 [infilename]			- .RData file containing phyloseq object from quick_and_dirty and make_phyloseq
											Default: paste("minimap_",primer_pair,"/Phyloseq_Outfile_MIMP_",primer_pair,".RData",sep="")
				3 [singleton_cutoff] 	- minimum number of sequences to carry forward from a given taxonomic assignment
				4 [taxon_level]			- the taxonomic level at which to aggregate and subset sequences 
				5 [rarefaction_level]	- the number of sequences to randomly subset from each taxonid
				6 [threads]				- the number of threads to run
				7 [kingdom_column]		- informs the script which column of the in file tax table contains kingdom information, for compatibility with samtools/minimap2 output from various databases
				8 [fasta_folder]		- folder to which to save the taxon-aggregated, rarefied fasta files
				9 [resume]				- a flag to tell the script to resume [?RESUME WHAT?] where it left off	
			Output: cluster_by_taxon/fasta_folder/cluster_xxx.fa

	Options:
		-p primer_pair='ITS54'		- the primer pair name (current default 'ITS1FLR3')
		-f psfilename=''			- .RData file containing phyloseq object from quick_and_dirty and make_phyloseq
		-l taxonlevel='genus'		- the taxonomic level at which to aggregate and subset sequences 
		-s singleton_cutoff=20		- minimum number of sequences to carry forward from a given taxonomic assignment
		-r rarefaction_level=1000	- the number of sequences to randomly subset from each taxonid
		-k kingdom_column=3			- informs the script which column of the in file tax table contains kingdom information, for compatibility with samtools/minimap2 output from various databases
		-t threads=8				- the number of threads to run
		-T tabulate_only='F'		- whether to only run cluster_by_taxon_p1_parallel.R, which produces a seqid-by-taxid table cluster_by_taxon/taxid_seqid.tsv
		-A aggregate_only='F'		- whether to only run cluster_by_taxon_p2_parallel.R, which does the aggregation by taxon and rarefaction
		-R fasta_folder='taxon_cluster' - folder to which to save the taxon-aggregated, rarefied fasta files
		-F resume='F'				- a flag to tell the script to resume [?RESUME WHAT?] where it left off	
		-L session_log_file="[-p primer_pair].log" 
				name of log file where you are keeping track of all the commands
				you are running and with which parameters (automatically supplies
				one based on primer pair)

	Output (via R scripts):
		cluster_by_taxon/taxid_seqid.tsv
		cluster_by_taxon/fasta_folder/cluster_xxx.fa
	
6. align_and_cluster_subtaxon.sh

	Usage: sh align_and_cluster_subtaxon.sh [-t subtaxon_alignment_threads] [-m mafft_threads] [-p phylip_threads] [-c mothur_threads] [-P mothur_path] [-A alignment only] [-D distance only] [-S subtaxon clustering only] [-B post subtaxon bin seqs only] [-C consensus sequences only] [-W pairwise alignment only] [-R resume task] [-i input directory to look for fasta files] [-o output directory for mothur] [-L session log file] [-h display this help message]

	Description: This script is a dynamic command that performs several tasks to ulimately cluster and generate consensus sequences:
		i. generates mafft alignment-based OR mothur pairwise.dist-based distance matrices
				i.a.1 mafft
				i.a.2 phylip dna.dist
					OR
				i.b mothur pairwise.dist
		ii. from distances in (i), "subtaxon" clustering with mothur cluster.classic [NOTE: need to implement directory specification]
		iii. "bin seqs" with mothur bin.seqs at the chosen distance cutoff label
		iv.  subsequent alignments of sequences within the subtaxon clusters and generation of consensus sequences
	
	Options (note some options are currently applied for different specifications in multiple steps):
	
		General Options:
			-P mothur_path='~/mothur/mothur'
			-i input_fasta_directory='taxon_cluster' 		- input directory for mothur to look for fasta files
			-o output_mothur_directory='distance_matrices' 	- directory where mothur saves output files (may need to change depending on steps being run)
			-A 												- only performs the first MAFFT alignment step (i.a.1)
			-D												- only performs the phylip dnadist step (i.a.2)
			-S												- only performs "subtaxon" clustering (ii) in mothur
			-B												- only runs mothur bin.seqs
			-C												- only does subtaxon alignments and consensus sequence calling
			-W												- generates pairwise distance matrices in mothur with pairwise.dist
			-L session_log_file=""							- specifies name of log file (default align_and_cluster_subtaxon.log)
	
		MAFFT and phylip alignment and distance options (i.a)
			-t subtaxon_alignment_threads=4	- number of instances of mafft to run
			-m mafft_threads=4				- number of threads to run in mafft
			-p phylip_threads=16			- number of instances of phylip to run
			-R resume='F'					- pick up where left off (currently only implemented here)
		
		mothur pairwise.dist Options
			-c mothur_threads=16			- number of instance of mothur to run
		
		Subtaxon clustering (mothur cluster.classic)
			-t subtaxon_alignment_threads=4	- number of instances of mothur to run
			[NOTE: need to implement directory specification]

	Output files:
		i.a.1	taxon_alignments/*.maf
		i.a.2	distance_matrices/*.maf.dist
		i.b		output_mothur_directory/*.phylip.dist
		ii.		*.rabund [NOTE: need to implement directory specification]
				*.sabund [NOTE: need to implement directory specification]
				*.list [NOTE: need to implement directory specification]
		iii.	*.an.[cutoff].fasta [NOTE: need to implement directory specification]
		iv.		*.cons.denovo.python.fasta	- final seq file
				*.cons.denovo.python.data	- data file
				*.reformatted.fasta			- seq file
				*.karuna.cons.ungap			- seq file
				*.karuna.cons				- seq file
				*.cons.data					- data file

c) MIMP "Sanger"
----------------


LAST UPDATED FEB 6 2024
