
import sys
import os
import subprocess as sp
import shutil
import numpy
from Bio import SeqIO

''''
Usage: python trim_filter.py [directory with all 'exon.fasta'  files]

sample 1 needs to be the outgroup
sample 2 needs to the better of the two ingroup transcripts. 
You will create a blastdb of sample1, and blast sample2 to sample1. 
------------------------
written for Python 2.7
Kyle O'Connell
kyle.oconnell@mavs.uta.edu
August 2016
------------------------
'''
#sample1 = 'maculata'
#sample2 = 'lacerataN'
#sample3 = 'lacerataS'


#define the base directory with the fasta files
fasta_directory = sys.argv[1]
#move to the base directory
os.chdir(fasta_directory)
#create directory to move the output files into
out_dir = fasta_directory+'/'+'Trimmed'
if not os.path.exists(out_dir):
	os.mkdir(out_dir)


sample1 = str(raw_input("Please type the name of sample 1: "))
sample2 = str(raw_input("Please type the name of sample 2: "))
sample3 = str(raw_input("Please type the name of sample 3: "))
##########################################################
#parse the original 3 seq alignment and write it to 3 new files
#for file in this directory
for filetype in os.listdir("."):
	if filetype.endswith('.fasta'):
		#open the file
		fh_fasta = open(filetype, 'rU')
		#use Bioseq to parse the fasta, then write the contents of each 
		#seq to a new file, named after the sample_ENSID
		for record in SeqIO.parse(fh_fasta, "fasta"):
			output_handle = record.id + '.fasta1'
			SeqIO.write(record, output_handle, "fasta")
		fh_fasta.close()
	
##########################################################
#Rewrite the files so that the seqs are in a single line
def fasta_text(x):
	#Create an empty dictionary that will be populated with 'names' as keys and 'sequences' as values
	IDseq_dict={}	
	for line in x:
		#Strip to remove all whitespace and newline characters
		line=line.strip()
		line = line.replace('-','')
		line = line.replace('!','')
		#Search for the '>' to signal sequence name
		if line.startswith('>'):
			#create variable equal to everything in the line except the '>'
			sample_name=line[1:]
			split_name = sample_name.split('_')
			sample = split_name[0]
			if sample.startswith(sample1):
				sample = sample[:len(sample1)]
			elif sample.startswith(sample2):
				sample = sample[:len(sample2)]
			EID = split_name[1]
			combined_name = sample + '_' + EID
			IDseq_dict[combined_name]=''
		else:
			#Add the subsequent sequence to the value (for each line not starting with '>')
			IDseq_dict[combined_name]+=line
	#####THIS IS THE FILE YOU WANT TO USE FOR BLASTN AND PARSING LATER#######	
	#####SAMPLE.FASTA1 CAN BE REMOVED#####	
	file_temp = '{}.fasta2'.format(combined_name)
	fh_temp = open(file_temp,'a')
    #go through our dictionary of names and sequences one entry at a time
	for seq in IDseq_dict:
		fh_temp.write('>'+seq+'\n'+IDseq_dict[seq]+'\n')
	fh_temp.close
	
#execute the above function on each file
for filetype4 in os.listdir('.'):
	#if the fasta is named accordingly:
    if filetype4.endswith('.fasta1'):
    	#create a temporary file handle to open it
        fas_temp=open(filetype4, 'r')
        #read the file lines into memory
        fas_temp_lines=fas_temp.readlines()
        #apply fasta_text function to this file
        fasta_text(fas_temp_lines)
        #once function is complete for this file, close it
        fas_temp.close()
        #loop continues for all appropriate remaining files
        
##########################################################
#define function to blast sample2 (query) to sample1 (db)
def blastn(x,y):
	db = sample1 + '_' + Ensemble_comb1 + '.fasta2'
	query = sample2 + '_' + Ensemble_comb1 + '.fasta2'
	outfile = Ensemble_comb1 + '.blastout'
	proc_blastn = sp.Popen(['blastn', '-db', db, '-query', query, '-out', outfile, '-outfmt', "6 qstart qend sstart send"])
	proc_blastn.wait()
		
##########################################################
#create a blast db for sample 1 for each exon	
#run the function above
for filetype2 in os.listdir("."):
	if filetype2.endswith('.fasta2'):
		#only worry about the new files
		names = filetype2.split(".")
		name = names[0]
		ID = name.split("_")
		Ensemble_comb1 = ID[1]
		if filetype2.startswith(sample1):
			proc_blastdb = sp.Popen(['makeblastdb', '-in', filetype2, '-dbtype', 'nucl'])
			proc_blastdb.wait()
		#run the blast command on sample 1 and 2. Sample 1 is db, sample 2 is query
		filetype2 = filetype2.split('.')
		name = filetype2[0]
		name = name.split('_')
		#call the blast function
		blastn(sample2,sample1)	
	

##########################################################
#parse blast output file 
for filetype3 in os.listdir("."):
	if filetype3.endswith(".blastout"):
		filename = filetype3.split('.')
		EnsembleID = filename[0]
		blastout = open(filetype3, 'r')
		#create temp dicts to write sequence to
		s1_dict = {}
		s2_dict = {}
		coord_list = []
		#define the ENSID
		#read line by line
		for line in blastout:
			#strip
			line = line.strip()
			#split file by tab
			line = line.split('\t')
			#define the sequence coordinates for file 1 and 2
			#query (sample2)
			coord_list.append(line[0])
			coord_list.append(line[1])
			coord_list.append(line[2])
			coord_list.append(line[3])
			#define the fasta2 files for sample 1 and 2
			file_s1 = sample1 + '_' + EnsembleID + '.fasta2'
			file_s2 = sample2 + '_' + EnsembleID + '.fasta2'
			#open both files within the loop
			fh_s1 = open(file_s1, 'r')
			fh_s2 = open(file_s2, 'r')
			for line in fh_s1:
				line = line.strip()
				if line.startswith('>'):
					#create variable equal to everything in the line except the '>'
					sample_name=line[1:]
					split_name = sample_name.split('_')
					sample = split_name[0]
					combined_name = sample + '_' + EnsembleID
					s1_dict[combined_name]=''
				else:
					#Add the subsequent sequence to the value (for each line not starting with '>')
					s1_dict[combined_name]+=line[int(coord_list[2]):int(coord_list[3])]
			for line in fh_s2:
				line = line.strip()
				if line.startswith('>'):
					#create variable equal to everything in the line except the '>'
					sample_name=line[1:]
					split_name = sample_name.split('_')
					sample = split_name[0]
					combined_name = sample + '_' + EnsembleID
					s1_dict[combined_name]=''
				else:
					#Add the subsequent sequence to the value (for each line not starting with '>')
					s1_dict[combined_name]+=line[int(coord_list[0]):int(coord_list[1])]	
		file_temp = EnsembleID + ".fasta3"
		fh_temp = open(file_temp, 'a')
		for seq in s1_dict:
			fh_temp.write('>' + seq + '\n' + s1_dict[seq[0:999]]+'\n')
		fh_temp.close()
##########################################################
#trim any sequence longer than 1000 with fastx_trimmer
for filetype4 in os.listdir('.'):
	if filetype4.endswith('.fasta3'):
		print "---------------------------------------------------"
		print "trimming sequences over 1000 bp from", filetype4, "using 'fastx_trimmer'", '\n'
		print "---------------------------------------------------"
		filetype_split = filetype4.split('.')
		file_prefix = filetype_split[0]+'.trimmed.fasta'
		program = 'fastx_trimmer'
		input_file = filetype4
		output = file_prefix
		proc_fastx_trimmer = sp.Popen([program, '-l 1000', '-i{}'.format(input_file), '-o{}'.format(output)])
		proc_fastx_trimmer.wait()	
#aln with muscle
		input_m = output
		file_prefix_m = filetype_split[0]+'.trimmed.aln.fasta'
		output_m = file_prefix_m
		proc_muscle = sp.Popen(['muscle', '-in', input_m, '-out', output_m])
		proc_muscle.wait()
##########################################################
#Rewrite the files so that the seqs are in a single line
def fasta_text(x):
	#Create an empty dictionary that will be populated with 'names' as keys and 'sequences' as values
	IDseq_dict={}	
	for line in x:
		#Strip to remove all whitespace and newline characters
		line=line.strip()
		#Search for the '>' to signal sequence name
		if line.startswith('>'):
			#create variable equal to everything in the line except the '>'
			sample_name=line[1:]
			split_name = sample_name.split('_')
			sample = split_name[0]
			EID = split_name[1]
			combined_name = sample + '_' + EID
			IDseq_dict[combined_name]=''
		else:
			#Add the subsequent sequence to the value (for each line not starting with '>')
			IDseq_dict[combined_name]+=line
	file_temp = '{}.final.fasta'.format(EID)
	fh_temp = open(file_temp,'a')
    #go through our dictionary of names and sequences one entry at a time
	for seq in IDseq_dict:
		fh_temp.write('>'+seq+'\n'+IDseq_dict[seq]+'\n')
	fh_temp.close
	
#execute the above function on each file
for filetype4 in os.listdir('.'):
	#if the fasta is named accordingly:
    if filetype4.endswith('.trimmed.aln.fasta'):
    	#create a temporary file handle to open it
        fas_temp=open(filetype4, 'r')
        #read the file lines into memory
        fas_temp_lines=fas_temp.readlines()
        #apply fasta_text function to this file
        fasta_text(fas_temp_lines)
        #once function is complete for this file, close it
        fas_temp.close()
        #loop continues for all appropriate remaining files
##########################################################

#clean up yo
for filetype5 in os.listdir('.'):
	if filetype5.endswith('.trimmed.fasta') or filetype5.endswith('.trimmed.aln.fasta') or filetype5.endswith('.fasta1') or filetype5.endswith('.fasta2') or filetype5.endswith('.fasta3') or filetype5.endswith('.nhr') or filetype5.endswith('nin') or filetype5.endswith('nsq') or filetype5.endswith('blastout'):
		os.remove(filetype5)
	if filetype5.endswith('.final.fasta'):
		proc = sp.Popen(['mv', filetype5, out_dir])
		proc.wait() 
#--------------------------------------------------------------------------------

print "-----------------------------------------------------------------------------------------------------"	   
print "Script is complete, please examine results in 'Fasta_Directory'", '\n'
print "-----------------------------------------------------------------------------------------------------"	
	

