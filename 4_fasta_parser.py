''''
Usage: python fasta_parser.py [directory with all processed fasta files 'Trimmed']

*Must first run gc_counter.py, and gc_count.out.txt must be in the dir with the trimmed files

This script will parse the gc_count info file, and then sort the processed fasta files accordingly and output two files:
1) 4-pruned_exon_list.txt 
2) 4-repeat_maker_upload.fasta

It will also print to the screen the number of kept and filtered exons. You can modify the parameters below change the rigor of filtering

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
import numpy as np

##############################################
#Modify to change filtering criteria

#user input for choosing low cutoff
lowcutoff = float('1')
#user input for choosing high cutoff
highcutoff = float('15')
gc_low =float('40')
gc_high=float('70')
length_low=float('400')
length_high=float('10000')

##############################################

#assign the directory where all the trimmed fasta files are
Trimmed_fasta_directory = sys.argv[1]
os.chdir(Trimmed_fasta_directory)

#open the gc info file and read it into memory
filename = 'gc_count.out.txt'
fh = open(filename, 'r')
#create output 1
filename2 = '4-pruned_exon_list.txt'
fh2 = open(filename2, 'a')
#write header
fh2.write('ENS_ID' + '\t' + 'Length' + '\t' + 'pdistance' + '\n')
#create file 2 and write header line
filename4 = '4-repeat_masker_upload.fasta'
fh4 = open(filename4, 'a')
fh4.write('#Repeat Masker Upload file: a concatenated file of all individuals and all loci' + '\n' + '\n')

#create a number of lists for future use
exon_list = []
pruned_exon_set = set()

#read the fh line by line to populate the lists
for line in fh:
	#skip the header
	if line.startswith('Ensemble'):
		pass
	else:
		#split by tab
		line_split = line.split('\t')
		#assign variables to each value in each column 
		data = line_split[0:4]
		exon_list.append(data)
#create empty dictionary
IDseq_dict = {}
#remember that entry = ID, gc, length, and pdist. It is the whole line from the gc info file
for entry in exon_list:
	Ens_ID = entry[0]
	gc = entry[1]
	length = entry[2]
	pdist = entry[3]
	#I am ditching the gc info, I only want to write the length and pdist to the filtered exon output
	combined_name = length + '\t' + pdist
	#filter the loci based on the parameters defined above
	if float(gc) > gc_low and float(gc) < gc_high and float(pdist) > lowcutoff and float(pdist) < highcutoff and float(length) > length_low and float(length) <= length_high:
		#if it passes, write the name and length/pdist to a dictionary
		IDseq_dict[Ens_ID]=combined_name

print "-----------------------------------------------------------"
print "total number of starting exons:", len(exon_list)
print
print 'total number of passed exons:', len(IDseq_dict)
print
print 'total number of exons that did not pass filtering:', len(exon_list) - len(IDseq_dict)
print "-----------------------------------------------------------"

#Create the Repeat Masker input file, if the filename equals the filtered stuff.
for filetype in os.listdir('.'):
	file_name = filetype.split('.')
	file_name = file_name[0]
	#match key to filename in dir to determine if the file is in the filtered dictionary keys
	for key in IDseq_dict:
		if key == file_name:
			#if the key and filename match, write it to a new file that includes ID length pdist
			fh2.write(key + '\t' + IDseq_dict[key] + '\n')
			#write all lines from the open fasta files to a new file, repeat maker input
			fh3 = open(filetype, 'r')
			for line in fh3:
				line = line.strip()
				line = line.replace('-','')
				fh4.write(line + '\n')
#close file for appending
fh4.close()			
#Find and print average length and pdist of retained loci based on filtering 
len_list = []
pdist_list = []

for value in IDseq_dict.values():
	#split the length and pdist within the value
	value = value.split('\t')
	length = value[0]
	pdist = value[1]
	#assign the length and pdist to differnt lists
	len_list.append(length)
	pdist_list.append(pdist)
	
#convert the lists to arrays to find the mean with numpy
len_array = np.asarray(len_list[1:]).astype(np.float)
pdist_array = np.asarray(pdist_list[1:]).astype(np.float)
#find the mean of each array
len_ave = np.mean(len_array)
pdist_ave = np.mean(pdist_array)

print "-----------------------------------------------------------"
print "The average length across all retained genes is {}bp".format(len_ave)
print "The average pdist across all retained genes is {}%".format(pdist_ave)
print "-----------------------------------------------------------"

#count the lines in the repeat masker input file
RM_lines = []
#open file for reading
fh4 = open(filename4, 'r')
#read line by line
for line in fh4:
	#strip of whitespace
	line = line.strip()
	#write the lines to a list
	RM_lines.append(line)
#divide the length by 2 to account for names and sequences
file_length = len(RM_lines)/2
#subtract one to remove the header line from the count
file_length = file_length - 1
	
#print the final length of the repeat masker file and the number of exons X 3. 
print "-----------------------------------------------------------"
print "Repeat Masker fasta upload written to {}. The file should have {} sequences".format(filename4, 3*len(IDseq_dict))
print
print "The actual number of sequences is {}".format(file_length)	
print
print "If these are different you may have an issue with the pipeline"
print
print "Run Repeat Masker on your output file, then move to script 5: Repeat_filter"
print "-----------------------------------------------------------"

fh.close()
fh2.close()	
fh3.close()
fh4.close()
