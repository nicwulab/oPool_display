
import glob
import multiprocessing

def is_quality_below_threshold(ascii_quality_scores, quality_threshold=40, max_low_quality_bases=5):
    quality_scores = [ord(char) - 33 for char in ascii_quality_scores]
    low_quality_bases = sum(qual < quality_threshold for qual in quality_scores)
    return low_quality_bases > max_low_quality_bases

def extract_and_filter_sequences(input_fastq_path, quality_threshold=40, max_low_quality_bases=5):
    split_filename = input_fastq_path.split("_")
    sample_name = split_filename[1][:-4]
    output_fastq_path = 'fastq_filtered/' + sample_name + '_filtered.fastq'
    print(output_fastq_path)
    with open(input_fastq_path, "r") as input_fastq, open(output_fastq_path, "w") as output_fastq:
        good_count = 0
        bad_quality_count = 0
        bad_length_count = 0
        not_found_count = 0
        while True:
            # Read four lines at a time (representing a single sequence entry in FASTQ)
            seq_id = input_fastq.readline().strip()
            sequence = input_fastq.readline().strip()
            qual_id = input_fastq.readline().strip()
            quality_scores = input_fastq.readline().strip()

            # Check for the end of file
            if not seq_id:
                break

            # Define multiple start/stop sequence pairs to search
            sequence_pairs = [
                ("GGAGTATCCACCATG", "TCCGGAGGATCCGAT"),
                ("ATCGGATCCTCCGGA", "CATGGTGGATACTCC"),
                # Add more sequence pairs as needed
            ]

            trimmed_sequence = None
            trimmed_quality_scores = None

            # Try each start/stop sequence pair for extraction
            for start_seq, stop_seq in sequence_pairs:
                start_index = sequence.find(start_seq)
                end_index = sequence.find(stop_seq) + len(stop_seq) if sequence.find(stop_seq) != -1 else -1

                # If either start or stop sequence is found within the read
                if start_index != -1 or end_index != -1:
                    # Extract the region
                    trimmed_sequence = sequence[start_index:end_index] if start_index != -1 else sequence
                    trimmed_quality_scores = quality_scores[start_index:end_index] if start_index != -1 else quality_scores
                    break  # Exit the loop if successful

            if trimmed_sequence is None or trimmed_quality_scores is None:
                not_found_count += 1
                #print("Desired start/stop sequences not found in read:", seq_id)
                continue

            if is_quality_below_threshold(trimmed_quality_scores, quality_threshold, max_low_quality_bases):
                bad_quality_count += 1
                #print("Low quality scores found in read:", seq_id)
                continue

            # Write to FASTQ file
            output_fastq.write(f"{seq_id}\n{trimmed_sequence}\n{qual_id}\n{trimmed_quality_scores}\n")
            good_count += 1

            #if good_count == 10:
                #break
    print(input_fastq_path)
    print('good count: ', good_count)
    print('bad quality count: ', bad_quality_count)
    print('not found count: ', not_found_count)


# Example usage:
#input_fastq_list = glob.glob('fastq/*.fastq')  # Replace with your FASTQ file path

input_fastq_list = ['fastq/Sample_35-Fwd_bc_01-Rev_bc_08.fastq',
                    'fastq/Sample_36-Fwd_bc_02-Rev_bc_08.fastq',
                    'fastq/Sample_37-Fwd_bc_03-Rev_bc_08.fastq',
                    'fastq/Sample_38-Fwd_bc_04-Rev_bc_08.fastq', 
                    'fastq/Sample_39-Fwd_bc_05-Rev_bc_09.fastq',
                    'fastq/Sample_40-Fwd_bc_06-Rev_bc_09.fastq']

# Replace with desired output FASTA file path
for input_fastq_path in input_fastq_list:
    extract_and_filter_sequences(input_fastq_path)
