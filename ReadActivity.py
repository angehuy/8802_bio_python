
'''
------------------------ Notes for myself ------------------------------------

- Run in terminal with: python ReadActivity.py Mzebra_GT3_Chr1.fasta
in the _GT3.fasta file, there are multiple accessin numbers which represent a different chromosome
So we need to iterate through all the values in the dictionary


- Need to add functionality for double RAD sequencing
'''

from pyfaidx import Fasta
import argparse
import re

parser = argparse.ArgumentParser(description = 'This script takes in a fasta file and an adaptor motif and identifies percentage of the fasta file that is sequenced') 
parser.add_argument('DNA_file', type = str, help = 'DNA sequence to analyze') # positional argument; storing it as a string type object
parser.add_argument('--Motif', type=str, help='Adaptor motif', default='GAATTC') # optional argument; needs -- or -


args = parser.parse_args()
file = args.DNA_file
motif = args.Motif

genes = Fasta(file) # creates a dictionary of accession:sequence
genes_names = genes.keys()
genes_names = list(genes_names) # gets a list of the accessions inside file

full_seq = {}


seq_adaptor_coord = []

total_genome_length = 0

for name in genes_names: #iterating through the ordered dict
    seq = genes[name] # gets the sequence
    seq = str(seq)
    seq_length = len(seq)
    total_genome_length += seq_length

    for match in re.finditer(motif, seq):
        coord = [match.start(), match.end()]
        seq_adaptor_coord.append(coord)

# checking if there are pairs of adaptor occurences; if not, then remove the last one
if len(seq_adaptor_coord) % 2 != 0:
    seq_adaptor_coord.pop()


# #1 = first adaptor
# #2 = second adaptor
# want to look at where #1 ends and #2 begins
# then for each of those, subtract the positions


inbetween_adaptor_count = 0


for i in range(0, len(seq_adaptor_coord), 2): # iterate through each pair of adaptors
    one_end = int(seq_adaptor_coord[i][1]) # getting #1's end
    two_start = int(seq_adaptor_coord[i+1][0])
    diff = two_start - one_end
    inbetween_adaptor_count += diff


percent_inbetween_adaptor = (inbetween_adaptor_count/ total_genome_length) * 100
percent_inbetween_adaptor_rounded = round(percent_inbetween_adaptor,2)

print(f'Percentage of genome sequenced by: {percent_inbetween_adaptor_rounded} %')





