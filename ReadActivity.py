'''
Example usage: python ReadActivity_class.py Mzebra_GT3_Chr1.fasta SingleRad

# This script currenlty has the capacity for SingleRad mode currently


The script works by:
1. Reading the FASTA file and creating a dictionary of accession IDs and their corresponding sequences.
2. Finding the locations of the adaptor motifs within each sequence.
3. Calculating the percentage of each sequence that is "sequenced" (based on the motifs), with different approaches depending on the mode (SingleRad or ddRad).
4. Returning the results as a dictionary, where the keys are accession IDs and the values are the calculated sequenced percentage for each accession.

Outputs:
- A dictionary mapping each accession ID to the percentage of the sequence that has been sequenced based on the provided adaptor motifs.
- A combined sequenced percentage across all accessions.

'''
from pyfaidx import Fasta
import argparse


parser = argparse.ArgumentParser(description = 'This script takes in a fasta file and an adaptor motif and identifies percentage of the fasta file that is sequenced') 
parser.add_argument('GenomeFile', type = str, help = 'DNA sequence to analyze') # positional argument; storing it as a string type object
parser.add_argument('Mode', type = str, choices = ['SingleRad', 'ddRad'], help = 'What type of RAD to run') # positional argument; storing it as a string type object
parser.add_argument('--re1', '--RE1', type=str, help='Restriction Enzyme / Adaptor motif', default='GAATTC')
parser.add_argument('--re2', '--RE2', type = str, help ='Optional argument to use with ddRAD mode') # optional argument; needs -- or -

args = parser.parse_args()
genome = args.GenomeFile
motif = args.re1


def readSequence(input_file: str):
    '''
    Parse a fasta file to create a accession:sequence dictionary
    '''
    seq = ''
    acc_seq = {}
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line[0] == ">": # If it's an accession line
                if seq:  # If sequence collected from previous run, collect it
                    acc_seq[id] = seq


                id = line[1:] # id stored in this line
                seq = ''  # Reset the sequence for the current accession
            else:
                seq += line  # Append the sequence data to the current sequence

        if seq: # For the very last line since there's no ">" to complete the GC calculations & dictionary adding
            acc_seq[id] = seq

        # print(type(acc_seq))
        return acc_seq


def find_motifs(seq: str, motif: str):
    '''
    Finds all the starting positions in a sequence and outputs into a dictionary of accession:(list of start positions)
    '''
    positions = []
    start = 0
    while True:
        start = seq.find(motif,start)
        if start == -1:
            break
        positions.append(start)
        start += len(motif)
    # print(positions)
    return positions


def calculate_seq_per_simple(seq: str, motif_positions: list):
    '''
    Calculates the percentage of sequence sequenced (simple assumption of 200 bases for each motif) 
    '''
    num_motifs = len(motif_positions)
    seq_length = len(seq)
    total_sequenced = num_motifs * 200
    # percentage_sequenced = (total_sequenced / seq_length) * 100
    return total_sequenced


def calculate_seq_per_extra(seq: str, motif_positions: list):
    '''
    Calculates the percentage of sequence sequenced (200 bases for each motif) 
    with 100bp flanking before and after the motif, if possible (catches edge cases if applicable)
    '''
    total_sequenced_bases = 0
    seq_length = len(seq)

    for start_position in motif_positions:
        # 100 bases before the motif and 100 after
        left_flank_start_pos = max(start_position - 100, 0) 
        # to ensure if position doesn't fully have 100bp ahead, if position - 100 results in a negative value 
        # (which happens when the motif is near the start of the sequence), the code will set the start position to 0 instead of a negative value
        right_flank_end_pos = min(start_position + 100 + len(motif), seq_length)
        # if position + 100 + len(motif) exceeds seq_length (the end of the sequence), it will be limited to seq_length
        left_flank = seq[left_flank_start_pos:start_position]
        right_flank = seq[start_position+len(motif) : right_flank_end_pos]
        sequenced_length = len(left_flank) + len(right_flank)
        
        total_sequenced += sequenced_length

    # Calculate percentage sequenced
    percentage_sequenced = (total_sequenced_bases / seq_length) * 100
    return total_sequenced


# Running for loop to create a giant dictionary of all the start positions of motif
acc_seq = readSequence(genome) # generates an accession:seq

motif_positions = {}
sequenced_percentages = {}
seq_len_total = 0
sequenced_len_total = 0

for acc,seq in acc_seq.items(): # for each acc:seq entry
    positions = find_motifs(seq,motif) # gets list of all positions
    motif_positions[acc] = positions # --> acc: [motif_pos1, motif_pos2]
    sequenced_percentages[acc] = calculate_seq_per_simple(seq, positions) # --> acc:length
    
    # combining across all accessions
    seq_len_total += len(seq)
    sequenced_len_total += sequenced_percentages[acc]

sequenced_percentages = round((sequenced_len_total / seq_len_total) * 100, 2)
print(f"The percentage of the genome sequenced (across all accessions) {sequenced_percentages}%.")














