''''
Usage: python 6.2-write_output.py [directory with all 'trimmed.fasta'  files]
Script will read in the final filtered exon output from script 5, and then write individual fasta files for each sample

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
fh1 = open(filename1, 'a')
fh2 = open(filename2, 'a')
#fh3 = open(filename3, 'a')

#create a list of ENSIDs from the filtered file output from script 5
ENSID_list = []

#open the filtered list of exons from script 5
filename4= '5-filtered_exon_list.txt'
fh4= open(filename4, 'r')

#read file line by line: ENSID, length, pdist
for line in fh4:
	#split line by tabs
	line_split = line.split('\t')
	#grab the 0 character, Ensemble ID
	ENS_ID = line_split[0]
	#write ENS_ID to ENSID_list
	ENSID_list.append(ENS_ID)
print
print
print "---------------------------------------------------------"
print "There are {} exons in the filtered exon list".format(len(ENSID_list))
print

#create empty dictionary where we will read in all sequences from the fasta ortholog files
IDseq_dict = {}

#iterate through our working dir. In this case 'Trimmed'
for filename in os.listdir('.'):
	#split the filename
	if filename.endswith('.final.fasta'):
		filename_split = filename.split('.')
		#Grab the Ensemble ID from the name.
		file_name = filename_split[0]
		#if Ensemble ID from filename is in the list of filtered exons from script 5
		if file_name in ENSID_list:
			#then open the file
			fh = open(filename, 'r')
			#read line by line
			for line in fh:
				#strip white spaces and newline char
				line = line.strip()
				#identify name lines
				if line.startswith('>'):
					#split by the carrot to grab sample name, including ENS ID
					line_split = line.split('>')
					#grab sample name
					name1 = line_split[1]
					#add sample name as key to dict
					IDseq_dict[name1]=''
				else:
					#next line is sequence
					seq = line
					#assign as dict value
					IDseq_dict[name1]+=seq

#Now we have a master dictionary of all sequences from filtered list. 
#create counts to determine the number of sequences in each file
i = 0 
j = 0
#k = 0
#iterate through our master dictionary
for key in IDseq_dict.keys():
	#determine if the samplename is in each key
	#if it's a hit, then count it, and write it to the appropriate file. 
	if key.startswith(str(sample1)):
		i = i + 1
		fh1.write('>'+key+'\n'+IDseq_dict[key]+'\n')
	if key.startswith(str(sample2)):
		j = j + 1
		fh2.write('>'+key+'\n'+IDseq_dict[key]+'\n')
	#if key.startswith(str(sample3)):
		#k = k + 1
		#fh3.write('>'+key+'\n'+IDseq_dict[key]+'\n')
	
print "The script is finished. Outputs were writen to '{}' and '{}'.'".format(filename1, filename2)
print
print "{} has {} exon sequences, and {} has {} exon sequences".format(filename1, i, filename2, j)
print "---------------------------------------------------------"
fh1.close()
fh2.close()
#fh3.close()