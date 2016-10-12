'''
Usage: python S2_choose_markers.py

*assumes 'S4_filtered_markers.txt' in working directory

Chooses a subset of markers from the 'S4_filtered_markers.txt' file
which should be in same directory.  There are three options for choosing a subset of
markers, listed below.  Program asks user input for choice, then
executes one of the three functions.  Output is a single file,
'S5_selected_markers.txt', containing ENSXETP # of marker, average
divergence, and length (with no headers).  File will be ranked by
divergence, lowest to highest.

*****************************
Option 1: "Random" - Generates a random distribution based on number of
markers requested, pulls out subset of markers contained in
distribution.  Input number of markers to choose. 

*****************************
Option 2: "Top" - Pulls a subset of genes with the highest divergence in file 
(skims off top divergence markers to use).  Input number of markers to choose.

*****************************
Option 3: "Combined" - Pulls a subset of genes with the highest
divergence in file, then from the remaining markers creates a random
distribution and pulls out subset of markers contained in the
distribution. Input number of top markers, then number of random
markers.  Combines both lists in output, ranked by divergence.

****Originally written by Dan Portik; Edited by Kyle O'Connell October 2016*****
'''

# **** Requires numpy v1.8 *****

import numpy as np

filename = '5-masked_exon_list.txt'
fh = open(filename, 'r')
fh2 = open('5-filtered_exon_list.txt','w')

#Create an enumerated dictionary list of the genes, divergences, and lengths (ie number all the genes from 1 to whatever in order of least to most divergent)

pruned_exon_list = []

#read through our file and get information we need from each line
lines = fh.readlines()
for line in lines:
    if (line.startswith('ENS') ):
        #will read in all lines beginning with "ENSXETP" (a gene)
        split_line=line.split()
        #Splits up each line into indices we use to put specific attributes into our list
        pruned_exon_list.append(split_line)
        #makes a list with three components: gene name, avgDiv, and Length

#intiate empy dictionary, we'll store the number of the marker as the key and all info from the list above as the value
exon_dict = {}

#assign numbers (keys) to the lists (values) (ie, [1:(ENSXETP, avgDiv, length)]) by enumerating and populating a dictionary
for i, x in enumerate(pruned_exon_list):
    exon_dict[i]=[]
    #add ENSXETP value to this key
    exon_dict[i].append(x[0])
    #add avgDiv value to this key
    exon_dict[i].append(x[1])
    #add length value to this key
    exon_dict[i].append(x[2])

#count number of genes to print at beginning of functions
gene_number = int(len(exon_dict))



#Option 1, "Random" Function - Generate a random distribution based on number of markers requested, pull out subset
#*****************************************************************
#define our function as option_1, let user know they've selected this
def option_1(x):
    print "Option 1 - Random selection"
    print
    print "Total number of data set: ", gene_number
    print
	
	#test to make sure a number is input for marker number
    marker_no1 = None
    while marker_no1 is None:
        try:
            marker_no1 = int(raw_input("How many random markers to select (ex. 200)? "))
        except ValueError:
            print("Try again")

    #Create random distribution of genes
    rando_unsorted_dist1 = np.random.choice(gene_number, marker_no1, replace=False) #0-High, size of output, without replacement
    #sort the random distribution from low value to high
    rando_dist1 = np.sort(rando_unsorted_dist1)
    print
    #show us the numbers in the distribution
    print "random distribution = ", rando_dist1
    print
    
    #For a number in the random distribution, find the key matching this number in our exon_dict dictionary
    random_exon_list1 = list()
    for k in rando_dist1:
    	#if it matches, put it in our new list above
        random_exon_list1.append(exon_dict[k])

    #Parse through our subset marker list for a bit more information to write the output file
    exon_list = []
    avgDiv_list = []
    Length_list = []
    for entry in random_exon_list1:
        gene_name = entry[0]
        avgDiv=entry[1]
        Length=entry[2]
        exon_list.append(gene_name)
        avgDiv_list.append(avgDiv)
        Length_list.append(Length)
        #write information to output file 'S2_selected_markers.txt'
        fh2.write(gene_name+'\t'+avgDiv+'\t'+Length+'\n')
        
    #a function that counts total number of base pairs in selected marker list
    def bp_counter(x):
        #start the counter at zero
        total_bases = int(0)
        #for every record of bp length in our Length_list, add the bp value to the counter 
        for i in x:
            total_bases += int(i)
        print "total number of base pairs = ", total_bases
        print
    bp_counter(Length_list)


#Option 2, "Top" Function - Pull a set of genes with the greatest divergence
#***************************************************************************
#define our function as option_2, let user know they've selected this
def option_2(x):
    print "Option 2 - Top fastest markers"
    print
    print "Total number of markers data set: ", gene_number
    print

    #test to make sure a number is input for marker number
    marker_no2 = None
    while marker_no2 is None:
        try:
            marker_no2 = int(raw_input("How many of the top diverging markers to select (ex. 200)? "))
        except ValueError:
            print("Try again")

    #Defines the top set of markers in a numerical distribution
    #first define the range of numbers by subtracting marker number desired from total number of markers
    min_range1 = int((gene_number)-(marker_no2))
    #create the range of numbers
    top_dist1 = np.arange(min_range1,(gene_number),1)
    #show us the range of numbers selected
    print "top distribution = ", top_dist1
    print

    #For a number in the top distribution, find the key matching this number in our exon_dict dictionary
    top_exon_list1 = list()
    for k in top_dist1:
    	#if it matches, put it in our new list above
        top_exon_list1.append(exon_dict[k])

    #Parse through our subset marker list for a bit more information to write the output file
    exon_list = []
    avgDiv_list = []
    Length_list = []
    for entry in top_exon_list1:
        gene_name = entry[0]
        avgDiv=entry[1]
        Length=entry[2]
        exon_list.append(gene_name)
        avgDiv_list.append(avgDiv)
        Length_list.append(Length)
        #write information to output file 'S2_selected_markers.txt'
        fh2.write(gene_name+'\t'+avgDiv+'\t'+Length+'\n')


    #counts total number of base pairs in selected marker list
    def bp_counter(x):
        #start the counter at zero
        total_bases = int(0)
        #for every record of bp length in our Length_list, add the bp value to the counter
        for i in x:
            total_bases += int(i)
        print "total number of base pairs = ", total_bases
        print
    bp_counter(Length_list)


#Option 3, "Combined" Function - Take a subset of markers with the greatest divergence, then random distribution of remaining markers
#***************************************************************************
#define our function as option_3, let user know they've selected this
def option_3(x):
    print "Option 3 - Subset of fastest markers, then random selection subset"
    print
    print "Total number of markers in data set: ", gene_number
    print

    #3a - Pull out top markers
    #first test to make sure a number is input for top marker number
    marker_no3 = None
    while marker_no3 is None:
        try:
            marker_no3 = int(raw_input("How many of the top diverging markers to select (ex. 200)? "))
        except ValueError:
            print("Try again")

    #Defines the top set of markers in a numerical distribution
    #first define the range of numbers by subtracting marker number desired from total number of markers
    min_range2 = int((gene_number)-(marker_no3))
    #create the range of numbers
    top_dist2 = np.arange(min_range2,(gene_number),1)
    print
    print "top distribution = ", top_dist2
    print

    #For a number in the top distribution, find the key matching this number in our exon_dict dictionary, write to list
    top_exon_list2 = list()
    for k in top_dist2:
        top_exon_list2.append(exon_dict[k])

    #3b - begin random selection from remaining markers
    #calculate remaining range of markers to choose from
    new_gene_number = int((gene_number)-(marker_no3))
    print "number of markers remaining = ", new_gene_number
    print

	#test to make sure a number is input for random marker number
    marker_no4 = None
    while marker_no4 is None:
        try:
            marker_no4 = int(raw_input("How many random markers to select from remaining markers (ex. 200)? "))
        except ValueError:
            print("Try again")

    #create random distribution of genes
    rando_unsorted_dist2 = np.random.choice(new_gene_number, marker_no4, replace=False) #0-High, size of output, without replacement
    rando_dist2 = np.sort(rando_unsorted_dist2)
    print
    print "random distribution = ", rando_dist2
    print

    #For a number in the random distribution, find the key matching this number in our exon_dict dictionary, write to list
    random_exon_list2 = list()
    for k in rando_dist2:
        random_exon_list2.append(exon_dict[k])

    #create final list of markers for Option 3 by combining the two lists
    combined_exon_list = random_exon_list2 + top_exon_list2

    #print marker number and list to screen
    print "total number of markers included = ", len(combined_exon_list)
    print

    #Parse through our subset marker list for a bit more information to write the output file
    exon_list = []
    avgDiv_list = []
    Length_list = []
    for entry in combined_exon_list:
        gene_name = entry[0]
        avgDiv=entry[1]
        Length=entry[2]
        exon_list.append(gene_name)
        avgDiv_list.append(avgDiv)
        Length_list.append(Length)
        #write information to output file 'selected_markers.txt'
        fh2.write(gene_name+'\t'+avgDiv+'\t'+Length+'\n')

    #counts total number of base pairs in selected marker list
    def bp_counter(x):
        total_bases = int(0)
        for i in x:
            total_bases += int(i)
        print "total number of base pairs = ", total_bases
    bp_counter(Length_list)

#******************************************************************************************************************


#User gets to define which method they would like to use for marker selection
user_choice = (raw_input("Select a method for marker selection (random, top, or combined): "))

#Carry out user choice by accessing above functions
if user_choice == "random":
    option_1(filename)
elif user_choice == "top":
    option_2(filename)
elif user_choice == "combined":
    option_3(filename)

print
print "All remaining filtered genes, avgDiv, and length printed to 'S5_selected_markers.txt'"
print
print "script complete"

#close our files when everything is written
fh.close()
fh2.close()
