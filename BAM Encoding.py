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
        print(record.id)
        print(record.seq)
        return record.seq[start:end] # After finding the right record, extract the correct portion of the DNA sequence depending on the start and end from the BET file

    return None

# Determine the keep rotating the state between the regions
def running_through_the_gap(current_state, start, end, bam_file):

    reads_in_the_gap = []

    with (pysam.AlignmentFile(bam_file, "rb") as bam):
        for read in bam:
            start_pos = read.pos  # This is 0-based in pysam
            querylength = 0
            for tuple in read.cigartuples:
                querylength += tuple[1]
            end_pos = querylength + start_pos - 1
            print("starting position of the read is", start_pos)
            print("ending position of the read is", end_pos)
            print("start of gap is", start)
            print("end of gap is", end)
            if ((start_pos) >= start and (start_pos) < end) or ((end_pos) >= start and (end_pos) < end) or ((start_pos) >= start and (end_pos) <= end) or ((start_pos) <= start and (end_pos) >= end):
                reads_in_the_gap.append(read)

        print("reads in the gap are:", reads_in_the_gap)

    gap_length = end - start
    for read in reads_in_the_gap:
        print("current read is", read)

        if read.cigartuples:
            print("the cigartuples are", read.cigartuples)
            inside_a_pileupcolumn = 0
            is_head = 0
            if current_state == 0: # State 0
                current_state = 1
                is_head = 1

            index_tuple = 0
            tuple_starting_position = read.pos
            cigartuples_rotating = {}

            for cigar in read.cigartuples:
                cigartuples_rotating[index_tuple] = tuple_starting_position
                tuple_starting_position += cigar[1]
                index_tuple += 1

            counting_tuple_index = 0

            for cigar in read.cigartuples:
                operation, length = cigar
                print(operation, length)

                if cigartuples_rotating[counting_tuple_index] + length - 1 < start:
                    counting_tuple_index += 1
                    continue
                elif cigartuples_rotating[counting_tuple_index] >= end:
                    counting_tuple_index += 0
                    continue
                elif cigartuples_rotating[counting_tuple_index] >= start:
                    position_in_tuple = 0
                    starting_tuple_position_record = cigartuples_rotating[counting_tuple_index]
                    counting_tuple_index += 1
                else:
                    position_in_tuple = start - cigartuples_rotating[counting_tuple_index]
                    starting_tuple_position_record = cigartuples_rotating[counting_tuple_index]
                    counting_tuple_index += 1

                print("current position in tuple is", position_in_tuple)

                if current_state == 1: # State 1
                    if operation == 0 or operation == 7 or operation == 8: # B
                        current_state = 1
                        counting_position = 0
                        # Rotating inside a cigar segment such as inside the 10M
                        while position_in_tuple <= length - 1:
                            if is_head == 1:
                                is_head = 0
                                position_in_tuple += 1
                                continue
                            if read.pos - start + position_in_tuple + inside_a_pileupcolumn >= gap_length: # End_of_a_region
                                break

                            counting_position += 1
                            position_in_tuple += 1

                    elif operation == 2 or operation == 3: # C
                        current_state = 2
                        counting_position = 0

                        print("inside_position", position_in_tuple)
                        print("inside_a_pileupcolumn", inside_a_pileupcolumn)
                        print("pileupcolumb.pos - start + inside_position + inside_a_pileupcolumn", read.pos - start + position_in_tuple + inside_a_pileupcolumn)

                        position_in_tuple += 1
                        inside_a_pileupcolumn += counting_position

                        current_state = 2

                        while position_in_tuple <= length - 1:
                            if read.pos - start + position_in_tuple + inside_a_pileupcolumn >= gap_length:
                                break
                            position_in_tuple += 1

                        inside_a_pileupcolumn += counting_position

                    elif operation not in [0, 7, 8, 2, 3]:
                        continue
                    else:
                        current_state = 3

                elif current_state == 2: # State 2
                    if operation == 2 or operation == 3: # E
                        current_state = 2
                        counting_position = 0
                        while position_in_tuple <= length - 1:
                            if read.pos - start + position_in_tuple + inside_a_pileupcolumn >= gap_length:
                                break

                            counting_position += 1
                            position_in_tuple += 1
                        inside_a_pileupcolumn += counting_position

                    elif operation == 0 or operation == 7 or operation == 8: # D
                        current_state = 1
                        counting_position = 0

                        counting_position += 1
                        position_in_tuple += 1
                        inside_a_pileupcolumn += counting_position

                        current_state = 1
                        counting_position = 0

                        # Rotating inside a cigar segment such as inside the 10M
                        while position_in_tuple <= length - 1:
                            if is_head == 1:
                                is_head = 0
                                position_in_tuple += 1
                                continue
                            if read.pos - start + position_in_tuple + inside_a_pileupcolumn >= gap_length: # End_of_a_region
                                break

                            counting_position += 1
                            position_in_tuple += 1

    return current_state


# Count transitions (A, B, C, D, E, F) from the BAM file, every strand
def count_transitions(bam_file, chrom, start, end, transition_counts, current_state, current_reads_in_the_region):
    # Initializing the transition_variable
    # Opening the bam_file
            region_length = end - start
            original_state = current_state
            for read in current_reads_in_the_region[chrom]: # The pileupread has two versions: the "forward version" & the "reverse version"
                """Forward and reverse strand counts"""
                current_state = original_state

                print("current read is", read)

                if read.is_reverse:
                    strand = 6 # Starting index for reverse transitions, first original, then reverse version
                elif read.is_forward:
                    strand = 0 # If the alignment is in forward, then the start index would be 0
                else:
                    continue
                print("the current strand is", strand)

                if read.cigartuples:
                    print("the cigartuples are", read.cigartuples)
                    inside_a_pileupcolumn = 0
                    is_head = 0
                    if current_state == 0: # State 0
                        current_state = 1
                        transition_counts[strand + 0][read.pos - start + inside_a_pileupcolumn] += 1
                        is_head = 1

                    index_tuple = 0
                    tuple_starting_position = read.pos
                    print("tuple starting position is", tuple_starting_position)
                    cigartuples_rotating = {}

                    for cigar in read.cigartuples:
                        cigartuples_rotating[index_tuple] = tuple_starting_position
                        tuple_starting_position += cigar[1]
                        index_tuple += 1

                    print("the cigar tuples dictionary is", cigartuples_rotating)
                    counting_tuple_index = 0

                    for cigar in read.cigartuples:
                        operation, length = cigar
                        print(operation, length)

                        if cigartuples_rotating[counting_tuple_index] + length - 1 < start:
                            print("cigar starting position", cigartuples_rotating[counting_tuple_index])
                            print("tuple length is", length)
                            print("start position is", start)

                            counting_tuple_index += 1
                            continue
                        elif cigartuples_rotating[counting_tuple_index] >= end:
                            counting_tuple_index += 0
                            continue
                        elif cigartuples_rotating[counting_tuple_index] >= start:
                            position_in_tuple = 0
                            starting_tuple_position_record = cigartuples_rotating[counting_tuple_index]
                            counting_tuple_index += 1
                        else:
                            position_in_tuple = start - cigartuples_rotating[counting_tuple_index]
                            starting_tuple_position_record = cigartuples_rotating[counting_tuple_index]
                            counting_tuple_index += 1

                        print("current position in tuple is", position_in_tuple)

                        if current_state == 1: # State 1
                            if operation == 0 or operation == 7 or operation == 8: # B
                                 current_state = 1
                                 counting_position = 0
                                 # Rotating inside a cigar segment such as inside the 10M
                                 while position_in_tuple <= length - 1:
                                     if is_head == 1:
                                         is_head = 0
                                         position_in_tuple += 1
                                         continue
                                     if read.pos - start + position_in_tuple + inside_a_pileupcolumn >= region_length: # End_of_a_region
                                         break
                                     transition_counts[strand + 1][starting_tuple_position_record + position_in_tuple - start] += 1
                                     counting_position += 1
                                     position_in_tuple += 1

                                 print("current counting matrix is \n", transition_counts)

                            elif operation == 2 or operation == 3: # C
                                 current_state = 2
                                 counting_position = 0

                                 print("inside_position", position_in_tuple)
                                 print("inside_a_pileupcolumn", inside_a_pileupcolumn)
                                 print("pileupcolumb.pos - start + inside_position + inside_a_pileupcolumn", read.pos - start + position_in_tuple + inside_a_pileupcolumn)
                                 transition_counts[strand + 2][starting_tuple_position_record + position_in_tuple - start] += 1
                                 position_in_tuple += 1
                                 inside_a_pileupcolumn += counting_position

                                 current_state = 2

                                 while position_in_tuple <= length - 1:
                                     if read.pos - start + position_in_tuple + inside_a_pileupcolumn >= region_length:
                                         break
                                 transition_counts[strand + 4][starting_tuple_position_record + position_in_tuple - start] += 1

                                 position_in_tuple += 1
                                 inside_a_pileupcolumn += counting_position

                                 print("current counting matrix is \n", transition_counts)


                            elif operation not in [0, 7, 8, 2, 3]:
                                continue
                            else:
                                current_state = 3
                                transition_counts[strand + 5][starting_tuple_position_record + position_in_tuple - start] += 1
                                print("current counting matrix is \n", transition_counts)


                        elif current_state == 2: # State 2
                            if operation == 2 or operation == 3: # E
                                current_state = 2
                                counting_position = 0
                                while position_in_tuple <= length - 1:
                                    if read.pos - start + position_in_tuple + inside_a_pileupcolumn >= region_length:
                                        break
                                    transition_counts[strand + 4][starting_tuple_position_record + position_in_tuple - start] += 1
                                    counting_position += 1
                                    position_in_tuple += 1
                                inside_a_pileupcolumn += counting_position
                                print("current counting matrix is \n", transition_counts)
                            elif operation == 0 or operation == 7 or operation == 8: # D
                                current_state = 1
                                counting_position = 0

                                transition_counts[strand + 3][starting_tuple_position_record + position_in_tuple - start] += 1

                                counting_position += 1
                                position_in_tuple += 1
                                inside_a_pileupcolumn += counting_position

                                current_state = 1
                                counting_position = 0

                                # Rotating inside a cigar segment such as inside the 10M
                                while position_in_tuple <= length - 1:
                                    if is_head == 1:
                                        is_head = 0
                                        position_in_tuple += 1
                                        continue
                                    if read.pos - start + position_in_tuple + inside_a_pileupcolumn >= region_length: # End_of_a_region
                                        break
                                    transition_counts[strand + 1][starting_tuple_position_record + position_in_tuple - start] += 1
                                    counting_position += 1
                                    position_in_tuple += 1

                                print("current counting matrix is \n", transition_counts)

            return transition_counts, current_state


def process_bed_fasta_bam(bed_file, fasta_file, bam_file, output_prefix):
    """Load the BED file using pd.read_csv"""
    bed_df = pd.read_csv(bed_file, sep='\s+', header = None, skiprows=1, names = ['chromosome', 'start', 'end', 'strand'])

    print(bed_df)
    last_state_of_region = {}

    current = 0

    current_reads_in_the_region = {}

    region_start = {}

    # breakpoint()
    for idx, row in bed_df.iterrows():
        print("the current row is", idx, row)

        chrom = row['chromosome']
        start = row['start']
        end = row['end']
        current_reads_in_the_region[chrom] = []
        region_start[chrom] = start

        with (pysam.AlignmentFile(bam_file, "rb") as bam):
            for read in bam:
                start_pos = read.pos  # This is 0-based in pysam
                querylength = 0
                for tuple in read.cigartuples:
                    querylength += tuple[1]
                end_pos = querylength + start_pos - 1
                print("starting position of the read is", start_pos)
                print("ending position of the read is", end_pos)
                print("start of region is", start)
                print("end of region is", end)
                if ((start_pos) >= start and (start_pos) < end) or ((end_pos) >= start and (end_pos) < end) or ((start_pos) >= start and (end_pos) <= end) or ((start_pos) <= start and (end_pos) >= end):
                   current_reads_in_the_region[chrom].append(read)
                print("current reads in the region " + chrom + " are", current_reads_in_the_region[chrom])
    # breakpoint()
    print(current_reads_in_the_region)

    sorted_dict = dict(sorted(region_start.items(), key=lambda item: item[1]))
    previous_ending = 0

    for idx, row in bed_df.iterrows():
        chrom = row['chromosome']
        start = row['start']
        end = row['end']

        print("the previous ending is", previous_ending)
        print("the start of the current region is", start)

        if (previous_ending < start):
            print("the ending state of the previous region is", current)
            current = running_through_the_gap(current, previous_ending, start, bam_file)
            print("now the correct start state of the current region is", current)

        ending_positon_of_region = end

        print("the start position is", start)
        print("the end position is, not included", end)

        length = end - start # the length of current region
        print("the length is", length)

        transition_counts = np.zeros((6, length)) # counting matrix for the current region

        sequence = extract_sequence_from_fasta(fasta_file, chrom, start, end)
        print("the start of the region is", start)
        print("the end of the region is", end)

        if sequence is None:
            continue

        # One-hot the sequence
        one_hot_encoded_sequence = one_hot_encoding(sequence)
        print("one hot encoded sequence", one_hot_encoded_sequence)

        # Get the 12 counts
        # breakpoint()
        transition_matrix, new_current = count_transitions(bam_file, chrom, start, end, transition_counts, current, current_reads_in_the_region)

        current = new_current
        last_state_of_region[chrom] = current

        # breakpoint()
        # Combine the one-hot encoding and the transition counts
        print("the transition counts matrix is \n", transition_matrix)
        print("the current one-hot encoding is", one_hot_encoded_sequence)

        combined_tensor = np.vstack([transition_counts, one_hot_encoded_sequence.T])
        print("the combined tensor is \n", combined_tensor)

        # Save the tensor as a numpy file
        output_file = f"{output_prefix}_region_{idx}.npy"
        np.save(output_file, combined_tensor)
        print(f"Save {output_file} with shape {combined_tensor.shape}")
        previous_ending = end

def main():
    # Define input files and output prefix
    bed_file = "sequence.bed"
    fasta_file = "sequence.fasta"
    bam_file = "sequence.bam"
    output_prefix = "output"
    sam_file = "sequence.sam"

    # Process the files
    process_bed_fasta_bam(bed_file, fasta_file, bam_file, output_prefix)

if __name__ == "__main__":
    main()
        

