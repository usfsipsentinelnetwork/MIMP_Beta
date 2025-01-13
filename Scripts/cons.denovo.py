import os
import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

file = sys.argv[1]

alignment = AlignIO.read(file, "clustal")
num_seqs = len(alignment)
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.gap_consensus(0.50, "-")

gapcons_data = str(consensus)

filename = file
gapcons_data = str(consensus)

with open("all_filt_denovo.maf.karuna.cons", "w") as gapcon:
    gapcon.write(">" + filename + "_gapped consensus" + "\n" + gapcons_data)

#print(summary_align)

with open("all_filt_denovo.maf.cons.data", "w") as gapcon_data:
    gapcon_data.write(str(summary_align.pos_specific_score_matrix()))
#This is where the gap characters get removed

with open ("all_filt_denovo.maf.karuna.cons", "r") as file:
	filedata = file.read()
	filedata_nogaps = filedata.replace("-", "") 
	filedata_ungap = filedata_nogaps.replace("_gapped ", "_")
 
with open("all_filt_denovo.maf.karuna.cons.ungap", "w") as file:
    file.write(filedata_ungap)
