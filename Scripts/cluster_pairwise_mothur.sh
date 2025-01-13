
~/mothur/mothur "#set.dir(output=/vol3/home/ec2-user/LM1-recalled-demultiplexed-80/cluster_by_taxon, input=/vol3/home/ec2-user/LM1-recalled-demultiplexed-80/cluster_by_taxon);pairwise.seqs(fasta=cluster_294.fa, output=lt, cutoff=100, processors=96)"

# dont need names file after all
# cd cluster_by_taxon
# Rscript -e "with(utils::read.table('taxid_seqid.tsv', header=1), utils::write.table(cbind(unique(seqid),unique(seqid)), 'seqids.names', quote=F, row.names=F, col.names=F))"

~/mothur/mothur "#set.dir(output=/vol3/home/ec2-user/LM1-recalled-demultiplexed-80/cluster_by_taxon, input=/vol3/home/ec2-user/LM1-recalled-demultiplexed-80/cluster_by_taxon);cluster.classic(phylip=cluster_294.phylip.dist, cutoff=100)"

#Rscript -e "with(utils::read.table('taxid_seqid.tsv' ..."

~/mothur/mothur "#set.dir(output=/vol3/home/ec2-user/LM1-recalled-demultiplexed-80/cluster_by_taxon, input=/vol3/home/ec2-user/LM1-recalled-demultiplexed-80/cluster_by_taxon);bin.seqs(list=cluster_294.phylip.an.list, fasta=cluster_294.fa, label=.22)"

ginsi --clustalout --inputorder --thread 64 --threadtb 20 --threadit 20 --quiet cluster_294.Otu1.fasta > cluster_294.Otu1.fasta.maf

python ../cons.denovo.py cluster_294.Otu1.fasta.maf

#../dnadist

#~/mothur/mothur "#set.dir(output=/vol3/home/ec2-user/LM1-recalled-demultiplexed-80/cluster_by_taxon, input=/vol3/home/ec2-user/LM1-recalled-demultiplexed-80/cluster_by_taxon);cluster.classic(phylip=cluster_294.fa.maf.dist)"

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < cluster_294.Otu1.fasta | sed '/^$/d' > cluster_294.Otu1.reformatted.fasta

Rscript ../consensus_seq_and_score.denovo.R cluster_294.Otu1.reformatted.fasta all_filt_denovo.maf.cons.data

##############
# 0. go back and filter out at higher cutoff for higher taxa (order 10000, class 100000 etc)
# 1. pairwise seq, output=lt
# 2. cluster.classic with cutoff = 100
# 3. needs to be some sort of analysis of cutoffs (can do in parallel with the genus-family level ones)
# 4. with chosen cutoff, bin.seqs at cutoff
# 5. separate sequences by bin
# 6. then align with MAFFT the binned OTUs!
# 7. call consensus sequences with python script
# 8. call consensus sequences and scores with mothur script
# 9. record all the consensus sequence scores
# 10.build a database of the consensus sequences
# 11.ref ids - use original UNITE classifications for original sequences and NCBI for finer scale taxonomy (for ~99-100% hits) and assign to reference sequences
# 12.rerun the alignment search (minimap) with the custom built database
# 13.quality check, check for low quality matches/hits - might need a second sweep/clean up

# need to work on documentation AND add more outputfiles (including a master output file and master config file) AND write up a protocol!

# default file name "pipeline-history.log to append command calls to"