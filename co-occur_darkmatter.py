##### This script is to figure out which dark matter (transcripts that had no diamond hits) co-occurs with the viral conserved segments
##### It takes four files at the command line: first a fasta file of the transcripts that did not have any blast hits, second a .clstr file from CD-hit analysis (cont.)
##### Third, a file that contains the transcript names of the viral conserved sequences. It should be a text file with one viral conserved seqence per line (no comma after each).
##### Fourth, it can take a diamond annotation file in .tsv format. This is to re-check if any of the transcripts have been annotated in other diamond runs. 
##### For this to work, the diamond file should be a tsv with the following fields (in this order): qseqid qlen sseqid qstart qend sstart send evalue bitscore length pident stitle

#Example call
#python3 co-occur_darkmatter.py dark_matter.fasta clusters_from_cdhit.clstr viral_transcripts_of_interest.txt diamond.annotations.dmnd.tsv 
#output is written to a file called co-occur_output.tsv. 

from Bio import SeqIO 
import os
import sys
from cdhit_reader import read_cdhit

#dark matter file
dm = list(SeqIO.parse(sys.argv[1], "fasta"))

#cd-hit cluster file
clstr = open(sys.argv[2], "r")

#initialize a dictionary for saving clusters
clust_dict = {}

#make clusters
clusters = read_cdhit(clstr).read_items()
#print(clusters)
for c in clusters:
	namelist = []
	for member in c.sequences:
		namelist.append(member.name.split("_")[0])
		clust_dict[c.name] = [*set(namelist)]
#print(clust_dict)

#initialize dark matter cluster dictionary
dm_clust = {}

#Get clusters from Dark matter
for tx in dm:
	dm_names = []
	for c in clusters:
		if tx.id == c.refname:
			for member in c.sequences:
				dm_names.append(member.name.split("_")[0])
				dm_clust[c.name] = [*set(dm_names)] 
			
		else:
			continue
	

#print(dm_clust)

#find clusters with viruses in it: This takes a list of the viral segments found in blast
toi_f = open(sys.argv[3], "r") #file with virus transcripts in it

#initialize list for transcript storage
toi = []

#make a list of the viral transcripts
for line in toi_f:	
	line = line.strip()
	toi.append(line)
#print(toi)

#initialize toi_dict
toi_dict = {}

#take the transcripts of interest and makes a dictionary of their cluster names and unique sample list
for line in toi:
	for c in clusters:
		toi_names = []
		if line == c.refname:
			for member in c.sequences:
				toi_names.append(member.name.split("_")[0])
				toi_dict[c.refname] = [*set(toi_names)]
#print(toi_dict)
#print(dm_clust)

co_occur = []

#create output for the screen
print("cluster" + "\t" + "Vco"+ "\t" + "Tco" + "\t" + "total_samps" + "\t" + "transcript"  ) #header
for k,v in toi_dict.items():
		print("\n")
		print("virus: " + k) 
		for l,w in dm_clust.items():
			cnt = 0
			for i in w:
				if i in v:
					cnt += 1	#the final cnt value is the intersection of V and T
				else:
					continue
			if cnt/len(v) >= 0.75 and cnt/len(set(w)) >= 0.5:
				  
				#print(l + "\t" + str(cnt/len(v)) + "\t" + str(len(w)/len(v)))
				for c in clusters:
					if  c.name == l:
						print(l + "\t" + str(cnt/len(v)) +  "\t" + str(cnt/len(set(w))) + "\t" + str(len(w)) + "\t" + c.refname) 
						co_occur.append([k,l, cnt/len(v),cnt/len(w),len(w), c.refname])
					else:
						continue	


#now check for annotations in dmnd file
dmnd_list = []
dmnd = open(sys.argv[4], "r")


for line in dmnd:
	line = line.strip()
	line = line.split("\t")
	dmnd_list.append(line)

for line in dmnd_list:
	for nl in co_occur:
		if line[-1] == nl[0]:
			nl.append(line[7])
			nl.append(line[11])
			break

#write out final output
of = open("co-occur_output.tsv", "w")
count = 0
for line in co_occur:
	if count == 0:	of.write("Viral Conserved Sequence" + "\t" + "Cluster" + "\t" + "Vco"+ "\t" + "Tco" + "\t" + "total_samps" + "\t" + "Transcript" + "Diamond Annotation" + "\n") #header

	if count > 0:
		of.write("\n")
	for element in line:
		of.write(str(element) + "\t")
	count += 1
	
		
#close files	
toi_f.close()
clstr.close()
dmnd.close()
of.close()


	



