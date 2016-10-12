'''
usage python 5-repeat_filter.py input_dir (Trimmed)

This script will cd into the Trimmed dir, grab the Repeat Modeler Output file, search for N's
then remove those loci from the list of final loci. 

*Requires that the repeat masker output file is the trimmed dir and is named 'repeat_masker_output.txt'
*Requires that pruned_exon_list.txt from fasta_parser.py is in the dir as well

------------------------
written for Python 2.7.3
Kyle O'Connell
kyle.oconnell@mavs.uta.edu
Sept 2016
------------------------

'''
import sys
import os
import subprocess as sp
import shutil
import numpy

#assign the directory where all the trimmed fasta files are
Trimmed_fasta_directory = sys.argv[1]
os.chdir(Trimmed_fasta_directory)

#open input and write the output file
filename = '4-repeat_masker_output.txt'
fh = open(filename, 'r')
filename2 = '4-pruned_exon_list.txt'
fh2 = open(filename2, 'r')
filename3 = '5-masked_exon_list.txt'
fh3 = open(filename3, 'a')

#make lists
mask_list = []
pruned_exon_list = []
final_list = []
#read through file line by line
for line in fh:
	#strip white space and newlines
	line = line.strip()
	#identify name
	if line.startswith('>'):
		line_split = line.split('_')
		Ens_ID = line_split[1]
	else:
		#identify seq lines
		seq = line
		#if N in seq searches for masked lines.
		if 'N' in seq:
			#if masked in any sequence, remove the whole exons, so all three species
			mask_list.append(Ens_ID)
# for line in pruned exon list
for line in fh2:
	#skip header line
	if line.startswith('ENS_ID'):
		pass
	else:
		#strip of newlines etc
		line = line.strip()
		#split by tab to sep. ID length and pdist
		line_split = line.split('\t')
		Ens_ID = line_split[0]
		#fill list of ENS IDs for comparison
		pruned_exon_list.append(Ens_ID)
		#Identify non masked sequences
		if Ens_ID not in mask_list:
			#write the non masked seqs to a new file
			fh3.write(line_split[0] + '\t' + line_split[1] + '\t' + line_split[2] + '\n')
			#fill final list to count the number of final exons easily
			final_list.append(Ens_ID)
difference = len(pruned_exon_list) - len(final_list)
		
print "----------------------------------------------------"
print "Number of masked sequences:", len(mask_list)
print
print "length of pruned list before masking:", len(pruned_exon_list)
print
print "length of final list after masking:", len(final_list)
print
print "{} exons were filtered out by masking".format(difference)
print
print "final list of {} filtered exons in output file '5-filtered_exon_list.txt'".format(len(final_list))
print "----------------------------------------------------"
fh.close()
fh2.close()
fh3.close()

