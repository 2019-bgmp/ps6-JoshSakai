#! /usr/bin/env python3
import argparse

def start_args():
    parser=argparse.ArgumentParser("Script for part 1 of PS6 to edit kmer files")
    parser.add_argument("-f", "--input_file", help="Input file", required=True)
    parser.add_argument("-k", "--kmer", help="holds kmer size", required=True, type=int)
    parser.add_argument("-o", "--output_file", help="output file, NOTE: use unique file names, output appends to files of the same name", required=True)
    return parser.parse_args()                  #A tragic lost opportunity for a pirate pun, args, arrrrrgs... man its 1am....
    
args = start_args()             #grabs the flags for the global statement
kmer = args.kmer                #holds the length of all the kmers
file = args.input_file                #sets the name of the input file in my naming convention of choice
output = args.output_file                #sets the name of the output file in my naming convention of choice

#for this PS kmer =49

#=====================================================================================================================================

with open(file, "r") as base:
    for line in base:
        line = line.strip("/n")         #BEGONE NEW LINE CHARACTER
        if line.startswith('>'):        #only the headers start with >
        #IMPORTANT NOTE, while working on testing this, delete old output files if you are using the same name, it will just add the outputs together.
        #hmm, should add that to the help in the argparse, done.
            with open(str(file)+"just_seq_id", "a") as build:        # I want to make a new file with just sequence ids
                build.write(line)
                
#=========================================================================================================================================
contig_length = []          #making an empty array to store kmers contig length
contig_coverage = []          #making an empty array to store kmers contig coverage

import re
with open(str(file)+"just_seq_id", "r") as working:
    for line in working:
        relevent = '(>)([A-Z]+)_([0-9]+)_([a-z]+)_([0-9]+)_([a-z]+)_([0-9, .]+)'                  #grabbing just the kmer_length of each contig and the kmer coverage for each contig
        element = re.search(relevent, line)                                                           #searches for what the regex specifies in the header
        if element:
            contig_length.append(float(element.group(5)))                       #code kept erroring out since it thought it was a string, made it a float
            contig_coverage.append(float(element.group(7)))                     #for the next part we only want contig_coverage


counter=0                   #initializing a counter to make the command iterative
adjust = kmer - 1                  #we want the actual contig length so we need to account for kmerizing
for element in contig_length:
    contig_length[counter]=element+adjust       #puts the actual contig lengths back into the list
    counter+=1

#==============================================================================================================================

#I need to define a lot of variables for this next part, it's killing me a little to do this peacemeal
#formulas taken from GenomeAssembly_and_Velvet lecture

#num_of_contigs is the number of reads in the input file
num_of_contigs = len(contig_length)

#does what it says on the tin, sugoi
max_contig_length = max(contig_length)

#thank you PS4
mean_contig_length = sum(contig_length)/num_of_contigs

#adding all the contigs lengths together gets you the total length of the genome 
genome_length = sum(contig_length)

#assuming all other calculations are correct, this should be the coverage.
#helpful advice from slides, if genome_length is unknown estimate from similar organism
#or, estimate it from kmer_freq in denovo genome projects 

kmer_coverage = []

con_count=0
for element in contig_length:
    coverage = (contig_coverage[con_count]*element)/(element-kmer+1)
    kmer_coverage.append(coverage)
    con_count+=1

mean_coverage = sum(kmer_coverage)/num_of_contigs


#================================================================================================================================
#calculating the N50

contig_length.sort(reverse= True)                   #sorting the values in the list from largest --> smallest
halfway_there = genome_length / 2                   #This tells me where the halfway point in my genome length is
#LIVING ON A P R A Y E R, rename variable something like pseudomedian for professional submission. OOOOEEEEE it's 2am.

N50_calc=0                      #initializing a counter to walk through the genome length until the halfway_there value
for element in contig_length:
    N50_calc+=element                   #counts through contig lengths in order of largest --> smallest
    if N50_calc < halfway_there:        #if the totals not the N50 value it moves right along
        pass
    elif N50_calc > halfway_there:           #breaks out of the loop as soon as the N50 is found, got the idea from PS5, thanks Leslie
        N50 = element                          #N50 = " at least half of the nucleotides in the assembly belongs to contigs with the N50 length or longer."
        break
        

#================================================================================================================================

buckets = {}            #for the final print statement we need to put the different contig lengths into buckets. I'm using a dictionary

for element in contig_length:
    ID = round(element, -2)         #rounds down to int
    if ID in buckets:
        buckets[ID]+=1
    elif ID not in buckets:
        buckets[ID]=1

#==================================================================================================================================

import operator  #Alright, I confess, I wasn't sure how to organize the dictionary things in a coherant print statement so I went to get a module... that's not too bad yeah? I did find it on my own and read the documentation
organize_buckets = sorted(buckets.items(), key=operator.itemgetter(0), reverse=False)                   #AHA, sorting the kmers into a tuple seems to have fixed my errors, I also needed to flip to sort from smallest --> largest
x=[]
y=[]
with open(output, "w") as out:
    print("Max Contig Length =", max_contig_length,file=out)
    print("Number of contigs =", num_of_contigs,file=out)
    print("Mean contig length =", mean_contig_length,file=out)
    print("Total genome length =", genome_length,file=out)
    print('N50 =', N50,file=out)
    print("Coverage =", mean_coverage,file=out)
    print("# Contig Length", "\t", "Number of Contigs in this Category",file=out)
    
    print_counter=0             #initializing a print counter to make this command iterative
    for element in organize_buckets:
        print(organize_buckets[print_counter][0], "\t", organize_buckets[print_counter][1],file=out)
        x.append(organize_buckets[print_counter][0])
        y.append(organize_buckets[print_counter][1])
        print_counter+=1
import matplotlib.pyplot as plt
plt.figure(figsize=(10,10))
plt.xlabel('contig length')
plt.ylabel('Number of contigs in this category')
plt.title("Distribution plot")
# plt.yscale('log')
# plt.xscale('log')
plt.scatter(x,y)
#print(x,y)
plt.savefig(str(file)+"_contig_distribution_"+str(kmer)+".png")
    #3am and FINISH

#back at it again at 1am, can't get this plot to show, I'm printing the values and it's fine. the axis are showing correct sizes
#after doing some googling it looks like the trouble is in the backend, going to install some packages and poke around

