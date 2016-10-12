''''
Usage: python gc_counter.py [directory with all 'trimmed.fasta' files]
*Best to run Aln_trim_transcripts.py first to remove gaps and trailing ends or poor alignements

This script will calculate the average GC, length, and pairwise distance for each locus alignment. Output is printed
to gc_count.out.txt and averages across all loci are printed to screen

*I commented out all references to sequence 3, since we trimmed our analysis down to two transcriptomes

------------------------
written for Python 2.7
Kyle O'Connell
kyle.oconnell@mavs.uta.edu
August 2016
------------------------
'''
#________________________________________#
import sys
import os
import subprocess as sp
import shutil
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Phylo.TreeConstruction import DistanceCalculator
import numpy
import numpy as np
#________________________________________#

#assign base dir and change to it, I recommend Trimmed
fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

#open the output file
filename = 'gc_count.out.txt'
fh = open(filename, 'a')
#write the header
fh.write('Ensemble_ID' + '\t' + 'Ave_GC' + '\t' + 'Ave_Length' + '\t' + 'Ave_pDist' + '\n')

#create lists
file_names = []
all_files = []
##################################################################
#print the number of files in the directory
#iterate through the trimmed dir and make a list of all file names
for filetype in os.listdir('.'):
    if filetype.endswith('.fasta'):
        all_files.append(filetype)
#count files in dir
file_no = len(all_files)

print
print "----------------------------------------------------------------------------"
print "There are {} fasta files in the current directory.".format(file_no)
print
##################################################################
#make a list of the total seq number, the gc and length of each sequence
#Open each file in base directory
for filetype in os.listdir('.'):
	#only grab fasta alignment files
    if filetype.startswith('ENS'):
    	fh2 = open(filetype, 'r')
        #split file on the '.' to split the name from trimmed.fasta
    	file_name_EID = filetype.split('.')
        #isolate the file name (ensemble ID)
        file_name = file_name_EID[0]
        #Append the EID to the list file_names. We will use this to write to the final output
        file_names.append(file_name)
        #create a series of new lists; some I don't use right now, all are within the loop so that I only write the three sequences from the file in question within the loop
        seq_list=[]
        gc_ave = []
        len_ave = []
        pdist_ave_list = []
        #read line by line for fasta file from dir
    	for line in fh2:
    		#Strip to remove all whitespace and newline characters
    		line=line.strip()
    		#Search for the '>' to signal sequence name
    		if line.startswith('>'):
    			#split by > to separate the > from the name
    			line_split = line.split('>')
    			name = line_split[1]
    			#name_list.append(name)
    		else:
                #if not > then the line is the sequence
    			seq = line
                #add the sequence to the list of sequences, which will reset each file, so it should never have more than 3 sequences in it
    			seq_list.append(seq)
                #read line by line in seq list
			for sequence in seq_list:
                #count gc % for each sequence
				gc_count = GC(sequence)
                #add GC values to a list (3 values each time it loops)
				gc_ave.append(gc_count)
                #calculate the length of each seq
				length = len(sequence)
                #append to length list
				len_ave.append(length)
	
##################################################################
		#now find the average between the two files for length and GC. Also find the pdistance
        #assign each sequence to a variable for the pdist calc.
        s1 = seq_list[0]
        s2 = seq_list[1]
        #s3 = seq_list[2]
        #calculate average of three gc %s
        GC_ave = sum(gc_ave)/len(gc_ave)
        #make it a string to write it to file
        GC_str_ave = str(GC_ave)
        #find average length (they are all the same)
        Len_ave = sum(len_ave)/len(len_ave)
        seq_len = Len_ave
        str_len = str(Len_ave)
        #===============================================#
        #this could be shortened with a function
        #calculate pairwise distance, and find average across the 3 seqs in each file, I multiple by 100 to make it a number rather than a decimal
        calculator = DistanceCalculator('identity')
        #pdist s1 and s2
        pd1 = (calculator._pairwise(s1,s2))*100
        #pdist s1 and s3
        #pd2 = (calculator._pairwise(s1,s3))*100
        #pdist s2 and s3
        #pd3 = (calculator._pairwise(s2,s3))*100
        #pd_ave = [pd1,pd2,pd3]
        #pd_mean = numpy.mean(pd_ave)
        #pd_mean2 = str(pd_mean)
        #write all outputs to the final output file
        fh.write(file_name + '\t' + GC_str_ave + '\t' + str_len + '\t' + str(pd1) + '\t'+'\n')
        fh2.close()
#close the file for appending, open below for reading     
fh.close()

#Part two-print summary statistics from the gc output file     
#*********************************************#        
#make lists for part 2
Ens_ID_list = [] 
gc_list = []
len_list = []
pdist_list = []

#open the file for reading
fh = open(filename, 'r')

#iterate through file line by line
for line in fh:
	if line.startswith('Ensemble'):
		pass
	else:
		#strip white space
		line = line.strip('')
		#split by tab
		line_split = line.split('\t')
		#assign variables to each value in each column 
		Ens_ID = line_split[0]
		gc = float(line_split[1])
		len = float(line_split[2])
		pdist = float(line_split[3])
		#append these values to the lists defined above
		Ens_ID_list.append(Ens_ID)
		gc_list.append(gc)
		len_list.append(len)
		pdist_list.append(pdist)
#conver the lists to arrays to find the mean with numpy
gc_array = np.asarray(gc_list[1:])
len_array = np.asarray(len_list[1:])
pdist_array = np.asarray(pdist_list[1:])

#find the mean of each array
gc_ave = np.mean(gc_array).astype(float)
len_ave = np.mean(len_array)
pdist_ave = np.mean(pdist_array)
max_length = max(len_list)
min_length = min(len_list)
max_pdist = max(pdist_list)
min_pdist = min(pdist_list)
print "The average GC content across all genes is {}%".format(str(gc_ave))

print "The average length across all genes is {}bp ({}-{}bp)".format(str(len_ave), str(min_length), str(max_length))

print "The average pdist across all genes is {}% ({}-{}%)".format(str(pdist_ave), str(min_pdist), str(max_pdist))
print "----------------------------------------------------------------------------"
print

fh.close()
