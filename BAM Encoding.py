import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO

# Define one-hot encoding for DNA, every base will be denoted by a specific vector of dimension equal to 4, respectively as 'A', 'T', 'G', 'C', 'N'
def one_hot_encoding(sequence):
    encoding = {
        'A': [1, 0, 0, 0],
        'T': [0, 1, 0, 0],
        'C': [0, 0, 1, 0],
        'G': [0, 0, 0, 1],
        'N': [0, 0, 0, 0]
    }
    return np.array([encoding[base] for base in sequence])

# Extract the DNA sequence based on the information from the BET file
def extract_sequence_from_fasta(fasta_file, chrom, start, end):
    for record in SeqIO.parse(fasta_file, "fasta"): # Using SeqIO to read through the record in the fasta file
        if record.id == chrom: # Using record.id and chrom to find the right record
            return record.seq[start:end] # After finding the right record, extract the correct portion of the DNA sequence depending on the start and end from the BET file
    return None

# Count transitions (A, B, C, D, E, F) from the BAM file, every strand
def count_transitions(bam_file, chrom, start, end, transition_counts, current_state):
    # Initializing the transition_variable
    # Opening the bam_file
    with pysam.AlignmentFile(bam_file, "rb") as bam: # Using pysam.AlignmentFile to open a file
        for pileupcolumn in bam.pileup(chrom, start, end): # Confirming the window Also using the start and end and chrom obtained from the BED file to extract relevant information from bam file
            # Every position(base) in the sequence will correspond to one row in the matrix
            for pileupread in pileupcolumn.pileups: # The pileupread has two versions: the "forward version" & the "reverse version"
                """Forward and reverse strand counts"""
                if pileupread.alignment.is_reverse:
                    strand = 6 # Starting index for reverse transitions, first original, then reverse version
                else:
                    strand = 0 # If the alignment is in forward, then the start index would be 0

                if pileupread.cigar:
                    for cigar in pileupread.cigar:
                        operation, length = cigar
                        if current_state == 0 and pileupread.is_head: # A
                            transition_counts[pileupcolumn.pos][strand + 0] += 1
                            current_state = 1
                        elif current_state == 1 and operation != 3: # B
                            transition_counts[pileupcolumn.pos][strand + 1] += 1
                            current_state = 1
                        elif current_state == 1 and operation == 3: # C
                            transition_counts[pileupcolumn.pos][strand + 2] += 1
                            current_state = 2
                        elif current_state == 2 and operation != 3: # D
                            transition_counts[pileupcolumn.pos][strand + 3] += 1
                            current_state = 1
                        elif current_state == 2 and operation == 3: # E
                            transition_counts[pileupcolumn.pos][strand + 4] += 1
                            current_state = 2
                        elif current_state == 1 and pileupread.is_tail: # F
                            transition_counts[pileupcolumn.pos][strand + 5] += 1
                            current_state = 3

    return transition_counts, current_state


def process_bed_fasta_bam(bed_file, fasta_file, bam_file, output_prefix):
    """Load the BED file using pd.read_csv"""
    bed_df = pd.read_csv(bed_file, sep='\t', header = None, names = ['chromosome', 'start', 'end'])

    current = 0
    for idx, row in bed_df.iterrows():
        chrom = row['chromosome']
        start = row['start']
        end = row['end']
        length = end - start + 1 # the length of current region

        transition_counts = [12][length] # counting matrix for the current region

        sequence = extract_sequence_from_fasta(fasta_file, chrom, start, end)

        if sequence is None:
            continue

        # One-hot the sequence
        one_hot_encoded_sequence = one_hot_encoding(sequence)

        # Get the 12 counts
        transition_matrix, new_current = count_transitions(bam_file, chrom, start, end, transition_counts, current)
        current = new_current

        # Combine the one-hot encoding and the transition counts
        combined_tensor = np.vstack([transition_counts, one_hot_encoded_sequence.T])

        # Save the tensor as a numpy file
        output_file = f"{output_prefix}_region_{idx}.npy"
        np.save(output_file, combined_tensor)
        print(f"Save {output_file} with shape {combined_tensor.shape}")
        

