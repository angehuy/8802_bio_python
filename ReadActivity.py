
'''
Given a restriction site, what percentage of the genome will be sequenced by
single RAD-seq?
 Given two restriction sites, what percentage of the genome will be sequenced
by double RAD-seq?
 Given a Type 2-RE site, what percentage of the genome will be sequenced?
'''

'''
Note: dont run during the play button just do in terminal

# Need the percentage of each base of your desired sequence
# Then multiply percentages based on the adaptor sequences


Usage: python ReadActivity.py Mzebra_GT3_Chr1.fasta

# want user to specify the motif

parser.add_argument('Motif', type = str, help = 'Adaptor motif')
motif = args.Motif

i.e. GAATTC

Are we trying to find the amount of bases within 2 GAATTC adaptor sequences?


---


in the _GT3.fasta file, there are multiple accessin numbers which represent a different chromosome
So we need to iterate through all the values in the dictionary
'''

from pyfaidx import Fasta
import argparse
import re

parser = argparse.ArgumentParser(description = 'This script takes in a fasta file') 
parser.add_argument('DNA_file', type = str, help = 'DNA sequence to analyze') # storing it as a string type object
args = parser.parse_args()
file = args.DNA_file

genes = Fasta(file) # creates a dictionary of accession:sequence
genes_names = genes.keys()
genes_names = list(genes_names) # gets a list of the accessions inside file
print(genes_names)

# print(genes['NC_036780.1'])

full_seq = {}

RAD_adaptor = "GAATTC"


for name in genes_names:
    seq = genes[name] # gets the sequence
    seq = str(seq)
    seq_length = len(seq)
    for match in re.finditer(RAD_adaptor, seq):
        print (match.start(), match.end())





# print(type(full_seq))


# print(genes.keys())
# Mzebra_Chr1 = genes['NC_036780.1']
# Mzebra_Chr1 = str(Mzebra_Chr1)
# print(type(Mzebra_Chr1))


# nuc_count = {
#     "A":0,
#     "T":0,
#     "G":0,
#     "C":0,
#     "N":0
# }



# for nucleotide in Mzebra_Chr1:
#     nuc_count[nucleotide] += 1

# print(nuc_count)


# G_per = (nuc_count["G"]) / len(Mzebra_Chr1)
# C_per = (nuc_count["C"]) / len(Mzebra_Chr1)
# A_per = (nuc_count["A"]) / len(Mzebra_Chr1 )
# T_per = (nuc_count["T"]) / len(Mzebra_Chr1)

# print(G_per)

# RAD_adaptor = "GAATTC"



# DB_adptor1 = "CGA"
# DB_adaptor2 = "TGC"


