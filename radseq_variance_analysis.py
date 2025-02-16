import argparse
import pyfaidx
import pdb
import vcf
from Bio.Restriction.Restriction_Dictionary import rest_dict

"""
Setup: 
- Create a conda environment for this homework and install pyfaidx, vcf, and biopython

Usage: 

- python radseq_variance_analysis.py Mzebra_GT3_Chr1.fasta YH_MC_samples_Chr1.vcf SingleRad EcoRI 
- python radseq_variance_analysis.py Mzebra_GT3_Chr1.fasta YH_MC_samples_Chr1.vcf ddRad EcoRI --RE2 HindIII


"""

def my_parse_args() -> argparse.Namespace:
    """
    Function uses argparse to parse the command-line options. This function has no arguments itself,
    but should return an argparse.Namespace object named "parsed_args" containing the command-line arguments.

    :return parsed_args: the parsed command-line arguments
    """

    parser = argparse.ArgumentParser(
        description='this script performs either SingleRad or ddRad sequencing on DNA from an input FASTA file and compared'
                    'the results to variants from an input VCF file')
    # TODO: Complete this function as described in the docstring
    parser.add_argument('GenomeFile', type = str, help = 'Path to the .fasta file containing the genome') 
    parser.add_argument('VCFFile', type = str, help = 'Path to the .vcf file with the DNA polymorphisms') 
    parser.add_argument('Mode', type = str, choices = ['SingleRad', 'ddRad'], help = 'What type of RAD to run') 
    parser.add_argument('RE1', type=str, help='Name of restriction enzyme')
    parser.add_argument('-re2','--RE2', type=str, help='Optional argment of name of second restriction enzyme to use in ddRad mode') # optional argument; needs -- or -


    parsed_args = parser.parse_args()
    return parsed_args


def read_fasta(path_to_fasta: str, chromosome: str) -> str:
    """
    Function reads the fasta file using pyfaidx. This function returns the dna sequence (as a str)
    associated with the specified chromosome.

    :param path_to_fasta: path to the fasta file, as a str
    :param chromosome: which chromosome to read. The default value of 'NC_036780.1' refers to the first chromosome
    :return dna: the dna associated with the specified chromosome, as a str
    """

    fasta_file = pyfaidx.Fasta(path_to_fasta) 
    dna = str(fasta_file[chromosome][:])  # Convert the sequence to string, without metadata
    return dna


def find_motifs(dna: str, motif: str): # -> list[int]
    """
    Function takes in DNA sequence (string) and a Restriction Enzyme motif (string) and will return a list of locations (integers) where the motif is found/starts

    :param dna: DNA sequence to be analyzed
    :param motif: restriction enzyme motif
    :return positions: list of ints of the motif positions
    """
    positions = []
    start = 0 
    while True:
        # pdb.set_trace()
        start = dna.find(motif,start) # .find() gets the first occurence of the motif; start -1 when there is no match
        if start == -1:
            break
        positions.append(start)
        start += len(motif)
    return positions


def run_single_rad(dna: str, re1: str): # -> list[tuple[int, int]]
    """
    Function takes in the dna (str) and the name of the restriction enzyme (str), 
    and returns a list of tuples where each tuple contains the start index (int) and stop index (int) of a
    site sequenced by SingleRad

    :param dna: DNA sequence to be analyzed
    :param re1: restriction enzyme name
    :return sequenced_sites: a list of tuples, where tuple contains the start and stop index of a sequenced site
    """

    seq_length = 100 # How long the sequencing reads are

    # Ensure the enzyme exists in the restriction dictionary
    assert re1 in rest_dict, f"Error: Restriction enzyme {re1} not found in rest_dict"

    # Get the recognition motif for the given enzyme
    motif = rest_dict[re1]['site']

    # Find the cut sites using the find_motifs function
    cut_sites = find_motifs(dna, motif)

    # Generate the sequenced sites with correct range
    sequenced_sites = [(site - seq_length, site + seq_length) for site in cut_sites]

    return sequenced_sites


def run_ddrad(dna: str, re1: str, re2: str): # -> list[tuple[int, int]]:
    """
    Function takes in the DNA sequence (str) and the names of the restriction enzymes, 
    and returns a list of tuples where each tuple contains the start index (int) and stop index (int)
    of a site sequenced by ddRad

    :param dna: DNA sequence to be analyzed
    :param re1: first restriction enzyme name
    :param re2: second restriction enzyme name
    :return sequenced_sites: a list of tuples, where each tuple contains the start and stop index of a sequenced site
    """
    min_size = 300  # Shortest distance for the DNA piece
    max_size = 700  # Longest distance for the DNA piece
    seq_length = 100  # How long the sequenced reads are

    # Verify that RE1 and RE2 are both in the rest_dict provided by biopython
    assert re1 in rest_dict, f"Error: Restriction enzyme {re1} not found in rest_dict"
    assert re2 in rest_dict, f"Error: Restriction enzyme {re2} not found in rest_dict"

    # Look up the motifs associated with the two enzymes
    motif1 = rest_dict[re1]['site']
    motif2 = rest_dict[re2]['site']

    # Use find_motifs function to find the sequencing sites for RE1 and RE2
    re1_sites = find_motifs(dna, motif1)
    re2_sites = find_motifs(dna, motif2)

    # Initialize an empty list to store sequenced sites
    sequenced_sites = []

    # Loop through each RE1 site
    for r1 in re1_sites:
        # Look for RE1 and RE2 sites left of the current RE1 site
        left_hits_r1 = [x for x in re1_sites if x < r1]
        left_hits_r2 = [x for x in re2_sites if x < r1]

        # Check conditions for left side
        if left_hits_r2 and min_size < (r1 - left_hits_r2[-1]) < max_size:
            if not left_hits_r1 or left_hits_r2[-1] > left_hits_r1[-1]:
                sequenced_sites.append((left_hits_r2[-1], r1))

        # Look for RE1 and RE2 sites right of the current RE1 site
        right_hits_r1 = [x for x in re1_sites if x > r1]
        right_hits_r2 = [x for x in re2_sites if x > r1]

        # Check conditions for right side
        if right_hits_r2 and min_size < (right_hits_r2[0] - r1) < max_size:
            if not right_hits_r1 or right_hits_r2[0] < right_hits_r1[0]:
                sequenced_sites.append((r1, right_hits_r2[0]))

    # To ensure we get both sides, check reverse of the logic too
    # Loop through each RE2 site
    for r2 in re2_sites: 
        # Look for RE1 and RE2 sites left of the current RE2 site
        left_hits_r1 = [x for x in re1_sites if x < r2]
        left_hits_r2 = [x for x in re2_sites if x < r2]

        # Check conditions for left side 
        if left_hits_r1 and min_size < (r2 - left_hits_r1[-1]) < max_size:
            if not left_hits_r2 or left_hits_r1[-1] > left_hits_r2[-1]:
                sequenced_sites.append((left_hits_r1[-1], r2))

        # Look for RE1 and RE2 sites right of the current RE2 site
        right_hits_r1 = [x for x in re1_sites if x > r2]
        right_hits_r2 = [x for x in re2_sites if x > r2]

        # Check conditions for right side
        if right_hits_r1 and min_size < (right_hits_r1[0] - r2) < max_size:
            if not right_hits_r2 or right_hits_r1[0] < right_hits_r2[0]:
                sequenced_sites.append((r2, right_hits_r1[0]))

    return sequenced_sites



def find_variable_sites(vcf_file_path: str, sequenced_sites):
    """
    Function takes in path to a VCF file and a list of sequenced sites (as returned by run_single_rad
    or run_ddrad) and return a list of sites that are both sequenced and contain variation

    :param vcf_file_path: Path of vcf file
    :sequenced_sites: List of tuples, where each tuple contains the start and stop index of a sequenced site, as returned 
    by your run_single_rad/run_ddrad()
    """

    # Read in the VCF file using the VCFReader class
    vcf_obj = vcf.Reader(filename=vcf_file_path)  # Open the VCF file

    sequenced_sites_variable = []  # List will contain all of the sequenced DNA pieces that have variation
    
    # Loop through every record
    for record in vcf_obj: 
        if record.CHROM != 'NC_036780.1':  # Only looking at records from Chromosome 1
            continue  # Skip non-chromosome 1 records
        
        # Ensure  both samples have called genotypes
        if record.num_called != 2:
            continue  # Skip records where both genotypes aren't called

        # Ensure one genotype is 0/0 and the other is 1/1 (check the genotype fields; use gt_nums?)
        genotypes = [sample['GT'] for sample in record.samples]
        if '0/0' in genotypes and '1/1' in genotypes:
            # If both if conditions pass, check if the variant's position falls within any sequenced site
            if any(start <= record.POS <= stop for start, stop in sequenced_sites): # if any of the vcf position in between start/stop of sequenced sites
                # If the variant falls within sequenced DNA, add the variant position and chromosome to the list
                sequenced_sites_variable.append((record.POS, record.CHROM))
    
    return sequenced_sites_variable




### MAIN CODE ###


def main():
    ### MAIN CODE ###
    args = my_parse_args()
    target_chromosome = 'NC_036780.1'
    dna_string = read_fasta(path_to_fasta=args.GenomeFile, chromosome=target_chromosome)

    if args.Mode == 'SingleRad':
        print(f'running SingleRad for chromosome {target_chromosome} and restriction enzyme {args.RE1}')
        seq_sites = run_single_rad(dna_string, args.RE1)
    elif args.Mode == 'ddRad':
        print(f'running ddRad for chromosome {target_chromosome} and restriction enzymes {args.RE1} and {args.RE2}')
        seq_sites = run_ddrad(dna_string, args.RE1, args.RE2)

    input_len = len(dna_string)
    n_sites = len(seq_sites)
    sequencing_length = seq_sites[0][1] - seq_sites[0][0]
    sequenced_fraction = sequencing_length * n_sites / input_len

    print(f'{args.Mode} sequencing complete.')
    print(f'\tlength of input sequence: {input_len}')
    print(f'\tnumber of sequencing sites located: {n_sites}')
    print(f'\tpercentage of nucleotides sequenced: {100 * sequenced_fraction:.3f}%')

    print('checking for variation within sequenced sites')
    variable_sites = find_variable_sites(args.VCFFile, sequenced_sites=seq_sites)
    n_var_sites = len(variable_sites)
    var_fraction = n_var_sites / n_sites

    print(f'\tnumber of sites with variation: {n_var_sites}')
    print(f'\tfraction of sites with variation: {var_fraction}')
    print('\nAnalysis Complete')


if __name__ == "__main__":
    main()






