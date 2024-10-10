import pandas as pd
import numpy as np
from Bio import SeqIO

def one_hot_encode_dna(sequence): # what codes will be used if
    """One-hot encode a DNA sequence. """
    encoding = {
        'A': [1, 0, 0, 0],
        'T': [0, 1, 0, 0],
        'C': [0, 0, 1, 0],
        'G': [0, 0, 0, 1],
        'N': [0, 0, 0, 0] # N represents any base
    }
    return np.array([encoding[base] for base in sequence])

def extract_sequecne_from_fasta(fasta_file, chrom, start, end):
    """Extract a sequence from the FASTA file based on the chromosome and position"""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == chrom:
            return str(record.seq[start:end])
    return None

def process_bed_and_fasta(bed_file, fasta_file, output_prefic):
    """From BED file get the chrom, start and end position,
    and then use these information to find the correponding DNA sequence
    and at last, one-hot encode the sequence"""
    """Load BED file, read the """
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chromosome", "start", "end"])

    for idx, row in bed_df.iterrows():
        chrom = row['chromesome']
        start = int(row['start'])
        end = int(row['end'])
        length = end - start + 1

        # Extratc the requested DNA seuqence using the information got from the previous operations and
        # then feed these information (chrome_name, start_pos, end_pos) to the function "extract_sequence_from_fasta"
        sequence = extract_sequecne_from_fasta(fasta_file, chrom, start, end)

        if sequence is None:
            print(f"Sequence not found for {chrom}:{start}-{end}")
            continue

        # one-hot encode the sequence
        one_hot_encode = one_hot_encode_dna(sequence)

        # Save the tensor as a numpy file
        output_file = f"{output_prefic}_region_{idx}.npy"
        np.save(output_file, one_hot_encode)
        print(f"Saved {output_file} with shape {one_hot_encode.shape}")

