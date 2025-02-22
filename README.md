# MIMP_Beta: MinIon Metabarcoding Pipeline

LAST UPDATED FEB 21 2025

```
++++++++++++++++++++++++++++++++++++++++++++
++                                        ++
++   MIMP: MinIon Metabarcoding Pipeline  ++
++   ===================================  ++
++                                        ++
++      (C) Geoffrey Williams 2025        ++
++                                        ++
++++++++++++++++++++++++++++++++++++++++++++
```

##### Recommended citation
>G. Williams, MIMP Beta: MinION Metabarcoding Pipeline (2025), GitHub repository, https://github.com/usfssentinelnetwork/MIMP_Beta

## I. About 

MIMP (MinIon Metabarcoding Pipeline) features executable shell scripts to which one can pass arguments to augment a basic metabarcoding pipeline for processing multiplexed amplicon data from ONT MinIon, after basecalling.

Currently the implementation is simple. MIMP only currently supports a basic high level (ie >= genus) amplicon pipeline

**1. "Quick & Dirty" which relies solely on alignment to the UNITE database.**

MIMP will also support:

**2. "De Novo" which uses preliminary UNITE alignments to inform within-taxon clustering**
(based on rarefied subsample of sequences) and re-alignment to resulting consensus sequences which can be BLASTed to a custom database (ie NCBI), for long-read amplicon metabarcoding. MIMP also supports
	> **NOTE: not yet fully implemented (bugs to resolve and final steps not yet complete):

MIMP is built from a preliminary version that can be used to assemble reads for Sanger-like purposes, while dealing with contaminants (which is not possible with Sanger). This functionality makes it well suited for sequencing from fungi in situ, provided a low level of contamination by other fungi:

**3. processing amplicon data from single-sample libraries:**
[MIMP "Sanger"](https://github.com/usfsipsentinelnetwork/MinIon_Sanger_Beta)

## II. Dependencies

For Workflow (1) - "quick and dirty"
------------------------------------
* Linux, a mac, or some platform where you can run bash (virtual machine, docker, etc.)
* [nanofilt](https://github.com/wdecoster/nanofilt) and [NanoPlot](https://github.com/wdecoster/NanoPlot)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](http://www.htslib.org/)
* [python](https://www.python.org/)   - check is likely already loaded/installed
* [R](https://www.r-project.org/) (and various R packages)
	* [dplyr](https://dplyr.tidyverse.org/) (NOTE: early versions needed dplyr-devel but should not any more)
	* [tidyr](https://tidyr.tidyverse.org/)


For Workflow (2) - "de novo" and Workflow (3): [MIMP "Sanger"](https://github.com/usfsipsentinelnetwork/MinIon_Sanger_Beta)
-----------------------------------

* [mafft](https://mafft.cbrc.jp/alignment/server/index.html)    - may need mpi version
* [mothur](https://mothur.org/)
* dnadist from [phylip](https://phylipweb.github.io/phylip/)
* Additional R packages
	* [Bioconductor](https://www.bioconductor.org/)
		* [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) (not on CRAN, BiocManager::install("BioStrings") - need BiocManager package)
		* [phyloseq](https://joey711.github.io/phyloseq/) (not on CRAN, BiocManager::install("phyloseq") - need BiocManager package)
	* [progress](https://github.com/r-lib/progress)
	* [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
	* [doSNOW](https://cran.r-project.org/web/packages/doSNOW/index.html)


Notes for MSU ICER users for installation on development nodes (2/20/2025)
--------------------------------------------------------------------------

You may run into issues loading the modules below but this should work on ICER (MSU). If not, check the output and try the different versions. This order worked on 2/20/2025 on devel node dev-amd24. Then you should be able to load the following (note these depend on gcc version 12.3). Run the following code in this order to load the dependencies
The first time you run this you should do this.
```
module load AOCC/4.0.0-GCCcore-12.3.0     #First, you are going to need to update the C compiler (for mafft)
module load Python/3.11.3-GCCcore-12.3.0  # and python
pip install NanoPlot
pip install nanofilt
pip install cutadapt
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2 # loads R, needed packages, and the BiocManager packages
module load minimap2/2.26-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0

#To add to the directory (careful, if you get this wrong you may have to log out and back in)
export PATH=$(pwd)/.local/bin:$PATH
```

In the future, you can put the rest in a script, without reinstalling the python packages. You'll need to load these every time you run the pipeline.
```
module load AOCC/4.0.0-GCCcore-12.3.0     #First, you are going to need to update the C compiler (for mafft)
module load Python/3.11.3-GCCcore-12.3.0  # and python
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2 # loads R, needed packages, and the BiocManager packages
module load minimap2/2.26-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0

#To add to the directory (careful, if you get this wrong you may have to log out and back in)
export PATH=$(pwd)/.local/bin:$PATH
```

>You'll also want to install the UNITE database or whatever reference database you need to supply to minimap2 (SILVA, custom refseq, etc). In the following instance, its in my home folder (~). 

That should be everything you need! If doing this on a job node, you will need to write a submission script. You can include the module load code in the script but may want to test it out first.

Needed only for workflows (2) - "de novo" & (3) "Sanger"
---------------------------------------------------
* All of the above
* For the MIMP "de novo" only (see below), you will need mafft and mothur. After the above, load them as follows...
	```
	module load MAFFT/7.520-GCC-12.3.0-with-extensions
	module load Mothur/1.48.0-foss-2023a-Python-3.11.3
	```
* Also needed: dnadist (from phylip - available as binary at [MinIon_Sanger_Beta](usfsipsentinelnetwork/MinIon_Sanger_Beta). You will need the executable, dnadist from phylip. Download it directly with the following command
	```
	wget -O dnadist https://github.com/usfsipsentinelnetwork/MinIon_Sanger_Beta/raw/refs/heads/main/Final_scripts/dnadist
	```

# III. Workflows 

## a) MIMP "Quick & Dirty"
-----------------------

### Description
MIMP "Quick & Dirty": output OTU table/ phyloseq object solely on alignment to the UNITE database. Suitable for aggregate genus level metabarcoding.

In order to run this version of the pipeline, you should download the repository and have the scripts in the parent directory containing your folders for each barcode. That is, you should have your output from basecalling arranged into folders "barcode01" or "barcode1" through "barcode**" and the following scripts/executables in the parent folder
```
primer_seqs.sh # with info on your primers
trim_and_sort.sh
quick_dirty_minimap.sh
get_minimap_output.sh -p ITS1F4
make_phyloseq.R ITS1F4 UNITE
dnadist
```

### Example protocol

1. First, you'll need to download this repository to your home directory (~) or wherever you are running the pipeline (HPCU, desktop, etc.). You can then download this directory
	```
	git clone https://github.com/usfsipsentinelnetwork/MIMP_Beta
	```

	> I recommend changing the name of the directory to your project\
	>
	>	```
	>	mv MIMP_Beta Your_project_directory_name
	>	```

2. Next, you'll want your sequences in individual folders within it. Assuming you are doing this from the same directory that that project is in,
	```
	cp Folder_containing_your_barcode_folders_with_fastq_files_in_them/barcode* Your_project_directory_name/Scripts/
	```

3. Then, you're going to want to make sure they are unzipped. If they aren't move into the directory and do it.
	```
	cd MIMP_Beta Your_project_directory_name/Scripts/
	```
	
	And then one of the following, depending on the extension. If they already end in fastq, you should be alright...
	```
	unzip barcode*/*.zip
	gzip -d barcode*/*.gz
	tar xvf barcode*/*.tar
	```

4. Make sure you edit primer_seqs.sh below if necessary, or find the name of the primer pair if its already in there to pass as an option (-p) to the bash scripts.

5. Run the pipeline. In this case, the primer sequences had been mostly trimmed off along with the adapters during basecalling, and there was very low quality, so we adjusted the options, filtering at Q10, very low minimimum sequence length for full ITS (-m 150), with no headcropping (-c 0) and only one base of "adapter overlap" (-o 1) which can incorrectly tag a sequence that didn't contain part of the primers used, but thats what we were dealing with in this particular case... For the most part I've only modified options that deviate from default settings. In the log file, you can see the entire call. Note that the Rscript step does not produce a log file.

	> NOTE: Typically, after the first step below (bash trim_and_sort.sh), you will want to look at the nanoplot results (.html files) and potentially rerun with adjusted options as needed. Also, as you can see from above, this pipeline might not produce results, and even throw errors at you, if you are using default settings with very low quality or incorrectly or wierdly preprocessed/basecalled/demuiltiplexed sequences. Always know what was done in terms of the kit used, basecalling methods, etc. to streamline troubleshooting.

	```
	bash trim_and_sort.sh -p ITS1F4 -q 10 -m 150 -c 0 -o 1
	bash quick_dirty_minimap.sh -p ITS1F4 -d ~/sh_general_release_dynamic_s_all_04.04.2024_dev.fasta
	bash get_minimap_output.sh -p ITS1F4
	Rscript make_phyloseq.R ITS1F4 UNITE
	```

6. Once you load this into R, you can generate a plot as follows:

	> NOTE: we are getting rid of columns in the taxonomy table that were introduced erroneously. Every database is different, and so the processing of minimap results can vary by database and database version. The script in this pipeline that reads parses the samtools output from minimap2 produces a phyloseq object where taxonomic ranks are named ta1, ta2, etc. We aggregated at the genus level, which is 'ta8' in this case...
	```
	library(dplyr)
	library(phyloseq)
	tax_table(psobj) <- tax_table(psobj)[,-c(1:2)]
	psobj %>%
		prune_taxa(rownames(tax_table(psobj))!='1',.) %>%
		microbiome::aggregate_taxa(level = 'ta8') %>%
		microViz::comp_barplot(tax_level = 'ta8')
	```

	This produced the following output
	![Example of figure from output](https://raw.githubusercontent.com/usfsipsentinelnetwork/MIMP_Beta/refs/heads/main/example_figure.png)

### 0. *primer_seqs.sh*
>This is not a command to run; this file contains all the information about primers and primer sequences that you may want to select from using the command line call; this script can be edited to add additional information in the same format that it is set up for GNU BASH langauge.

#### Usage:
NA

#### Description:

Each primer set (F and R) is given a number; here is the format used to specify seq, rev comp seq, and length min and max parameters:
```
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
```

### 1. *trim_and_sort.sh*

#### Usage:

```
bash trim_and_sort.sh [-p primer_pair] [-q quality_cutoff] [-m min_length] [-M max_length] [-c head_crop] [-a adapter_error] [-o adapter_overlap] [-N nanoplot_all(T/F)] [-n nanoplot_each(T/F)] [-x cpu_cores] [-s skip_to_cut_adapt] [-L session_log_file]
```

#### Description:
This script does the following
1. runs nanoplot (summary and visualization of quality and length, distribution etc. after basecalling)
	> NOTE: Typically, after trim_and_sort.sh, you will want to look at the nanoplot results (.html files) and potentially rerun with adjusted options as needed. If you are using default settings with very low quality or incorrectly or wierdly preprocessed/basecalled/demuiltiplexed sequences, you may get errors. Always know what was done in terms of the kit used, basecalling methods, etc. to streamline troubleshooting. At that point you can use the -x option to save time and skip nanoplot and nanofilt to head straight to cutadapt. But thats assuming you don't need to adjust nanofilt options, which include head crop (-c), quality filtering (-q), and min and max length (-m and -M). These are all parameters that viewing the results from nanoplot can help you choose.

2. runs nanofilt (filters by q score, sequence length, and trims some reads like adapters from the reads)
3. converts file to fasta format and shortens names of sequences in fasta sample headers for downstream compatibility
4. makes a new folder for the results for the primer set being used
5. runs cutadapt to trim forward and reverse primer seqs and reorients them and contatenates all seqs from each barcode

#### Options (and default values)

##### General
```
-p primer_pair='ITS1F4'    - the primer pair you want to use
-x skip_to_cut_adapt='F'    - whether to skip nanoplot & nanofilt and go straight to cutadapt
-L session_log_file="\[-p primer_pair\].log" 
	name of log file where you are keeping track of all the commands
	you are running and with which parameters (automatically supplies
	one based on primer pair)
```

##### NanoFilt Options
```
-q quality_cutoff=15    - the quality cutoff for nanofilt
-m min_length=200       - the minimum length for nanofilt
-M max_length=2000      - the maximum read length for nanofilt
-c head_crop=50         - the number of leading bases for nanofilt to trim
```

##### CutAdapt Options
```
-e adapter_error=.1    - adapter error rate
-o adapter_overlap=15  - the amount of overlap allowed
-x cpu_cores=8         - number of cores to use for cut adapt
```

##### NanoPlot Options
```
-N nanoplot_all='T'  - whether to make a nanoplot for all seqs combined
-n nanoplot_each='T' - whether to make a nanoplot for each barcode
```

##### Output files:
```
barcodeXX/nanoplot files
barcodeXX/all_filt_concatenated.fasta (before cutadapt)
barcodeXX/primer_pair/all_filt_reorient.fasta
```

### 2. *quick_dirty_minimap.sh*

#### Usage:
```
bash quick_dirty_minimap.sh [-p primer_pair] [-q quality_cutoff] [-d database] [-i infile] [-b barcode_threads] [-t minimap_threads] [-m minlen] [-s skip_minimap] [-S skip_samtools]  [-L session_log_file]
```

#### Description:
This script does the following
1. runs minimap2 (a mapping alingment search) for each barcode to the database of choice and produces samtools output
2. creates the samtools output summary for each barcode
3. runs process_minimap2.R to convert that into a OTU table\
	the R script parses the samtools output to make it readable in R,\
	particularly in make_phyloseq.R

#### Options (* required) (= default)

##### General Options
```
-p (=primer_pair)    - name of primer pair as it appears in primer_seqs.sh (Default 'ITS54')
-q (=quality_cutoff) - the quality cutoff (as phred score) to filter initial sequences ((Default 10, corresponds to 90%)
-m minlen=300        - minimum length of seqs to filter initial sequences
-s skip_minimap      - (Default 'F') can be used to skip the minimap step in combination with -S to skip the samtools summary of minimap step - note that skipping minimap but not samtools makes little sense unless you wanted to just do one step at a time; but the summary section will still run at the end to gather minimap output.
-S skip_samtools     - (Default 'F') can be used in combination with -s to skip minimap and processing of minimap results. effectively, this function only gathers and summarizes minimap results.
-L session_log_file="[-p primer_pair].log" 
	name of log file where you are keeping track of all the commands
	you are running and with which parameters (automatically supplies
	one based on primer pair)
```

##### MiniMap2 Options
```
-d* (=database)                      - the name and location of reference database file for taxonomic assignments (formatted as a fasta file)
-i infile='all_filt_reorient.fasta'  - name of input fasta file in each barcode folder
-b barcode_threads=4                 - number of minimaps to run independently (total cpus = b * t)
-t minimap_threads=4                 - number of minimap threads to tell minimap run in each minimap instance
-P minimap_path='minimap2'
```

#### Output files:
```
barcodeXX/primer_pair/all_filt_minimap.sam (minimap2 output)
barcodeXX/primer_pair/all_filt.samview.tsv (samtools view output)
barcodeXX/primer_pair/relabund_phred_q10.csv (database accession - like an OTU - raw abundance table)
```

### 3. *get_minimap_output.sh*

#### Usage:
```
bash get_minimap_output.sh [-p primer_pair] [-L session_log_file]
```
	
#### Description:
This reads in the minimap csv output from each folder and copies them all to a named file in the folder minimap_primer_pair
	
##### Options (* required)
```
-p* (=primer_pair) - name of primer pair as it appears in primer_seqs.sh
-L session_log_file="[-p primer_pair].log" 
	name of log file where you are keeping track of all the commands
	you are running and with which parameters (automatically supplies
	one based on primer pair)
```

##### Output:
```
./minimap_primer_pair/*.barcodeXX.csv
```

### 4. *make_phyloseq.R*
>** required for future steps, but does not write to session log file

#### Usage:
```
Rscript make_phyloseq.R [primer_pair*] [db*] [outfilename] [silva_path]
```

#### Description:
Uses preliminary quick and dirty minimap2 to create an OTU and taxon table to make a phyloseq object

#### Options (* required)
```
1 [primer_pair] - the primer pair name
2 [db]          - the database (currently supports =UNITE and =SILVA)
3 [outfilename] - name of .RData output file with phyloseq object
4 [silva_path]  - ***required if db=SILVA
```

#### Output:
```
outfilename.RData
```

## b) MIMP "De Novo"
--------------------
>**NOTE: not yet fully implemented (bugs to resolve and final steps not yet complete)

MIMP "De Novo": use preliminary UNITE alignments to inform within-taxon clustering (based on rarefied subsample of sequences) and re-alignment to resulting consensus sequences which can be BLASTed to a custom database (ie NCBI), for long-read amplicon metabarcoding

### 5. *aggregate_minimap_bytaxon.sh*
>**NOTE: not yet fully implemented (bugs to resolve and final steps not yet complete)

#### Usage:

```
bash aggregate_minimap_bytaxon.sh [-p primer_pair] [-f phyloseq filename] [-l taxon level for clustering] [-h print this help message] [-s singleton_cutoff] [-r rarefaction_level] [-k kingdom_column] [-t threads] [-T tabulate only] [-A aggregate only] [-R resume] [-F fasta_folder]
```

#### Description:
Uses output from preliminary quick and dirty minimap2 to do run the following R scripts\
(supplying the options given to it to the R scripts); descriptions and specs for R scripts below\

#### Options:
```
-p primer_pair='ITS54'      - the primer pair name (current default 'ITS1FLR3')
-f psfilename=''            - .RData file containing phyloseq object from quick_and_dirty and make_phyloseq
-l taxonlevel='genus'       - the taxonomic level at which to aggregate and subset sequences 
-s singleton_cutoff=20      - minimum number of sequences to carry forward from a given taxonomic assignment
-r rarefaction_level=1000   - the number of sequences to randomly subset from each taxonid
-k kingdom_column=3         - informs the script which column of the in file tax table contains kingdom information, for compatibility with samtools/minimap2 output from various databases
-t threads=8                - the number of threads to run
-T tabulate_only='F'        - whether to only run cluster_by_taxon_p1_parallel.R, which produces a seqid-by-taxid table cluster_by_taxon/taxid_seqid.tsv
-A aggregate_only='F'       - whether to only run cluster_by_taxon_p2_parallel.R, which does the aggregation by taxon and rarefaction
-R fasta_folder='taxon_cluster' - folder to which to save the taxon-aggregated, rarefied fasta files
-F resume='F'                   - a flag to tell the script to resume [?RESUME WHAT?] where it left off	
-L session_log_file="[-p primer_pair].log" 
		name of log file where you are keeping track of all the commands
		you are running and with which parameters (automatically supplies
		one based on primer pair)
```

#### Output (via R scripts below):
```
cluster_by_taxon/taxid_seqid.tsv
cluster_by_taxon/fasta_folder/cluster_xxx.fa
```

#### 5a. *cluster_by_taxon_p1_parallel.R*
>**NOTE: not yet fully implemented (bugs to resolve and final steps not yet complete)

##### Usage:
```
Rscript cluster_by_taxon_p1.R [primer_pair] [infilename] [kingdom_column] [threads]	
```

##### Description:
Tabulates seqids by taxon and outputs them to a file

##### Options:
```
1 [primer_pair]        - the primer pair name (current default 'ITS1FLR3')
2 [infilename]         - .RData file containing phyloseq object from quick_and_dirty and make_phyloseq
                          Default: paste("minimap_",primer_pair,"/Phyloseq_Outfile_MIMP_",primer_pair,".RData",sep="")
3 [kingdom_column]     - informs the script which column of the in file tax table contains kingdom information, for compatibility with samtools/minimap2 output from various databases
4 [threads]            - number of threads to run
```

##### Output:
```
cluster_by_taxon/taxid_seqid.tsv
```
			
#### 5b. *cluster_by_taxon_p2_parallel.R*
>**NOTE: not yet fully implemented (bugs to resolve and final steps not yet complete)

##### Usage:
```
Rscript cluster_by_taxon_p2.R [primer_pair] [infilename] [singleton_cutoff] [taxon_level] [rarefaction_level] [threads] [kingdom_column] [fasta_folder] [resume]
```

##### Description:
Uses output from p2 to aggregate sequences at the desired taxonomic and rarefaction levels

##### Options:
```
1 [primer_pair]         - the primer pair name (current default 'ITS1FLR3')
2 [infilename]          - .RData file containing phyloseq object from quick_and_dirty and make_phyloseq
                            Default: paste("minimap_",primer_pair,"/Phyloseq_Outfile_MIMP_",primer_pair,".RData",sep="")
3 [singleton_cutoff]    - minimum number of sequences to carry forward from a given taxonomic assignment
4 [taxon_level]         - the taxonomic level at which to aggregate and subset sequences 
5 [rarefaction_level]   - the number of sequences to randomly subset from each taxonid
6 [threads]             - the number of threads to run
7 [kingdom_column]      - informs the script which column of the in file tax table contains kingdom information, for compatibility with samtools/minimap2 output from various databases
8 [fasta_folder]        - folder to which to save the taxon-aggregated, rarefied fasta files
9 [resume]              - a flag to tell the script to resume [?RESUME WHAT?] where it left off	
```

##### Output:
```
cluster_by_taxon/fasta_folder/cluster_xxx.fa
```

### 6. *align_and_cluster_subtaxon.sh*
>**NOTE: not yet fully implemented (bugs to resolve and final steps not yet complete)

#### Usage:
```
bash align_and_cluster_subtaxon.sh [-t subtaxon_alignment_threads] [-m mafft_threads] [-p phylip_threads] [-c mothur_threads] [-P mothur_path] [-A alignment only] [-D distance only] [-S subtaxon clustering only] [-B post subtaxon bin seqs only] [-C consensus sequences only] [-W pairwise alignment only] [-R resume task] [-i input directory to look for fasta files] [-o output directory for mothur] [-L session log file] [-h display this help message]
```

#### Description:
This script is a dynamic command that performs several tasks to ulimately cluster and generate consensus sequences:
1. generates mafft alignment-based OR mothur pairwise.dist-based distance matrices
	1 .a.1 mafft
	1 .a.2 phylip dna.dist
		OR
	1 .b mothur pairwise.dist
2. from distances in (i), "subtaxon" clustering with mothur cluster.classic [NOTE: need to implement directory specification]
3. "bin seqs" with mothur bin.seqs at the chosen distance cutoff label
4.  subsequent alignments of sequences within the subtaxon clusters and generation of consensus sequences
			
#### Options
>(note some options are currently applied for different specifications in multiple steps):
			
##### General Options:
```
-P mothur_path='~/mothur/mothur'
-i input_fasta_directory='taxon_cluster'
     - input directory for mothur to look for fasta files
-o output_mothur_directory='distance_matrices'
     - directory where mothur saves output files (may need to change depending on steps being run)
-A   - only performs the first MAFFT alignment step (i.a.1)
-D   - only performs the phylip dnadist step (i.a.2)
-S   - only performs "subtaxon" clustering (ii) in mothur
-B   - only runs mothur bin.seqs
-C   - only does subtaxon alignments and consensus sequence calling
-W   - generates pairwise distance matrices in mothur with pairwise.dist
-L session_log_file=""
     - specifies name of log file (default align_and_cluster_subtaxon.log)
```		

##### MAFFT and phylip alignment and distance options (i.a)
```
-t subtaxon_alignment_threads=4 - number of instances of mafft to run
-m mafft_threads=4              - number of threads to run in mafft
-p phylip_threads=16            - number of instances of phylip to run
-R resume='F'                   - pick up where left off (currently only implemented here)
```

##### mothur pairwise.dist Options
```
-c mothur_threads=16  - number of instance of mothur to run
```
				
##### Subtaxon clustering (mothur cluster.classic)
```
-t subtaxon_alignment_threads=4	- number of instances of mothur to run
	[NOTE: need to implement directory specification]
```

#### Output files:
```
i.a.1   taxon_alignments/*.maf
i.a.2   distance_matrices/*.maf.dist
i.b     output_mothur_directory/*.phylip.dist
ii.     *.rabund [NOTE: need to implement directory specification]
        *.sabund [NOTE: need to implement directory specification]
        *.list [NOTE: need to implement directory specification]
iii.    *.an.[cutoff].fasta [NOTE: need to implement directory specification]
iv.     *.cons.denovo.python.fasta  - final seq file
        *.cons.denovo.python.data   - data file
        *.reformatted.fasta         - seq file
        *.karuna.cons.ungap         - seq file
        *.karuna.cons               - seq file
        *.cons.data                 - data file
```

## c) [MIMP "Sanger"](https://github.com/usfsipsentinelnetwork/MinIon_Sanger_Beta)
-------------------

A pipeline for assembly of sanger sequences from MinIon. The original pipeline is published [here](https://github.com/usfsipsentinelnetwork/MinIon_Sanger_Beta)

LAST UPDATED FEB 21 2025
