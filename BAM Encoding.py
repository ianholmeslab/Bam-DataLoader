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
def count_transitions(bam_file, chrom, start, end):
    """ Initialize the count matrix for forward and reverse"""
    # Opening the bam_file
    with pysam.AlignmentFile(bam_file, "rb") as bam: # Using pysam.AlignmentFile to open a file
        for pileupcolumn in bam.pileup(chrom, start, end): # Also using the start and end and chrom obtained from the BED file to ectract relevant information from bam file
            # Every position(base) in the sequence will correspond to one row in the matrix
            pos = pileupcolumn.pos - start # relevant position with the requested column
            if 0 <= pos < (end-start):
                for pileupread in pileupcolumn.pileups: # The pileupread has two versions: the "forward version" & the "reverse version"
                    """Forward and reverse strand counts"""
                    if pileupread.alignment.is_reverse:
                        transition_type = 6 # Starting index for reverse transitions, first original, then reverse version
                    else:
                        transition_type = 0 # If the alignment is not in reverse, then the start index would be 0

                    # Count various transitions and update the count matrix accordingly ?????? How to decode the transition conditions??
                    # And for every version of the strand, there will be a checking to fill its corresponding row in the counting matrix depending the following checking
                    if pileupread.is_del:
                        transition_counts[pos][transition_type + 0] += 1
                    elif pileupread.is_refskip:
                        transition_counts[pos][transition_type + 1] += 1
                    elif pileupread.indel > 0:
                        transition_counts[pos][transition_type + 2] += 1
                    else:
                        transition_counts[pos][transition_type + 3] += 1
                        if pileupread.alignment.query_sequence[pileupread.query_position] != pileupread.reference_base:
                            transition_counts[pos][transition_type + 4] += 1 # Mismatch
                        transition_counts[pos][transition_type + 5] += 1

    return transition_counts


def process_bed_fasta_bam(bed_file, fasta_file, bam_file, output_prefix):
    """Load the BED file using pd.read_csv"""
    bed_df = pd.read_csv(bed_file, sep='\t', header = None, names = ['chromosome', 'start', 'end'])

    for idx, row in bed_df.iterrows():
        chrom = row['chromosome']
        start = row['start']
        end = row['end']
        length = end - start + 1

        sequence = extract_sequence_from_fasta(fasta_file, chrom, start, end)

        if sequence is None:
            continue

        # One-hot the sequence
        one_hot_encoded_sequence = one_hot_encoding(sequence)

        # Get the 12 counts
        transition_counts = count_transitions(bam_file, chrom, start, end)

        # Combine the one-hot encoding and the transition counts
        combined_tensor = np.hstack([one_hot_encoded_sequence, transition_counts])

        # Save the tensor as a numpy file
        output_file = f"{output_prefix}_region_{idx}.npy"
        np.save(output_file, combined_tensor)
        print(f"Save {output_file} with shape {combined_tensor.shape}")
        

