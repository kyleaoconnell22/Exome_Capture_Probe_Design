''''
Usage: python 7-validation.py [directory with all 'trimmed.fasta'  files]
Script will read in the final .fasta files for each sample from script 6, write the exons to a set, and then check to make sure there are not 
unique exons in any one set.

------------------------
written for Python 2.7
Kyle O'Connell
kyle.oconnell@mavs.uta.edu
August 2016
------------------------
'''
import sys
import os
import subprocess as sp
import shutil
import numpy

#define the base directory with the fasta files
fasta_directory = sys.argv[1]
#move to the base directory
os.chdir(fasta_directory)

#Define the sample names, in this case there are three, if more than three, major script modification is needed
sample1 = str(raw_input("Please type the name of sample 1: "))
sample2 = str(raw_input("Please type the name of sample 2: "))
#sample3 = str(raw_input("Please type the name of sample 3: "))

#Define the output filenames, samplename.fasta, I identify them as string, though no necessary here
filename1 = str(sample1) + '.fasta'
filename2 = str(sample2) + '.fasta'
#filename3 = str(sample3) + '.fasta'

#open the output files
fh1 = open(filename1, 'r')
fh2 = open(filename2, 'r')
#fh3 = open(filename3, 'r')

#define sets for each file's Ensemble IDs
set1 = set()
set2 = set()
#set3 = set()

#read file line by line to take ENS ID from file 1
for line in fh1:
	#strip white space and newlines
	line = line.strip()
	#identify name
	if line.startswith('>'):
		line_split = line.split('_')
		Ens_ID = line_split[1]
		set1.add(Ens_ID)
	else:
		pass
		
#same loop for file2
for line in fh2:
	#strip white space and newlines
	line = line.strip()
	#identify name
	if line.startswith('>'):
		line_split = line.split('_')
		Ens_ID = line_split[1]
		set2.add(Ens_ID)
	else:
		pass
'''
#same loop for file 3
for line in fh3:
	#strip white space and newlines
	line = line.strip()
	#identify name
	if line.startswith('>'):
		line_split = line.split('_')
		Ens_ID = line_split[1]
		set3.add(Ens_ID)
	else:
		pass

'''	
#begin our analysis of the "sets" we created
#combine all sets, which will not include duplicates (property of the "set"), and print how many unique 'ENSXET' names are present
#********** if adding additional files, add " | fh#_IDs_set " to the end of this next line
combined_set = set1 | set2
print
print "combined set of all unique markers = ", len(combined_set)
print

#subtract only 'ENSXET' names shared across every set, leaving any unshared names present, then print how many unshared 'ENSXET' names are present
#********** if adding additional files, add " ^ fh#_IDs_set " to the end of this next line
differences = set1 - set2
print "number of differences found across sets = ", len(differences)
print

#decide whether files pass the test, that all markers here must be present in every fasta file with no remainders
if len(differences) == 0:
    print "All markers across files are matched, good job!"
    print
else:
    print "WARNING: there are some markers not shared across all fasta files"
    print



fh1.close()
fh2.close()
#fh3.close()