from collections import Counter
from rich.progress import track
from Bio import AlignIO
import pandas as pd
import numpy as np
import itertools
import random
import time, os
import multiprocessing
import argparse
import sys

# codon from: https://www.biologicscorp.com/tools/CodonUsage

C_Table = {"UUU":"F","UCU":"S","UAU":"Y","UGU":"C",
"UUC":"F","UCC":"S","UAC":"Y",
"UUA":"L","UCA":"S","UAA":"*","UGA":"*",
"UCG":"S","UGG":"W",
"CUU":"L","CCU":"P","CAU":"H","CGU":"R",
"CCC":"P","CGC":"R",
"CCA":"P","CAA":"Q",
"CCG":"P","CAG":"Q",
"AUU":"I","ACU":"T","AAU":"N","AGU":"S",
"AUC":"I","ACC":"T","AAC":"N","AGC":"S",
"ACA":"T","AAA":"K",
"AUG":"M","ACG":"T","AAG":"K",
"GUU":"V","GCU":"A","GAU":"D","GGU":"G",
"GUC":"V","GCC":"A","GAC":"D","GGC":"G",
"GUA":"V","GCA":"A","GAA":"E","GGA":"G",
"GUG":"V","GCG":"A","GAG":"E","GGG":"G"}

def read_F2D(Fasta_DB):
    # Sequence file from csv format to dictionary format
    with open(Fasta_DB, 'r') as F:
        Tmp = F.readlines()
    Seq_dict = {}
    for i in range(len(Tmp)):
        if Tmp[i].startswith('>'):
            Seq_id = Tmp[i][1:].strip().split(' ')[0] # the first one is ID
            Seq_dict[Seq_id] = Tmp[i+1].strip()
    return Seq_dict

def Dict_update(Seq_dict, TB):
    # update the DNA seq from the table
    Seq_dict_update = Seq_dict.copy()
    
    # Debug: Print CSV columns and first few rows
    print(f"CSV columns: {TB.columns.tolist()}")
    print(f"CSV has {len(TB)} rows")
    if len(TB) > 0:
        print(f"First few Name values: {TB['Name'].head().tolist() if 'Name' in TB.columns else 'Name column not found'}")
    
    successful_updates = 0
    failed_updates = 0
    
    for id in Seq_dict.keys():
        id_short = id.split(':')[0]
        try:
            TB_subset = TB[TB.Name == id_short]
            if len(TB_subset) == 0:
                print(f"Warning: No match found for ID '{id_short}' in CSV")
                failed_updates += 1
                continue
                
            # Check if required columns exist
            required_cols = ['VH_AA', 'VL_AA', 'CDRH1_AA', 'CDRH3_AA', 'CDRL3_AA', 'VH_nuc', 'VL_nuc']
            missing_cols = [col for col in required_cols if col not in TB.columns]
            if missing_cols:
                print(f"Error: Missing required columns in CSV: {missing_cols}")
                failed_updates += 1
                continue
            
            VH_AA = TB_subset.VH_AA.tolist()[0]
            VL_AA = TB_subset.VL_AA.tolist()[0]
            CDRH1_AA = TB_subset.CDRH1_AA.tolist()[0]
            CDRH3_AA = TB_subset.CDRH3_AA.tolist()[0]
            CDRL3_AA = TB_subset.CDRL3_AA.tolist()[0]
            VH_nuc = TB_subset.VH_nuc.tolist()[0]
            VL_nuc = TB_subset.VL_nuc.tolist()[0]
            AA_combine = VL_AA + " " + VH_AA
            Nr_combine = VL_nuc + " " + VH_nuc
            Seq_dict_update[id] = {"AA": AA_combine, "Nr": Nr_combine, "CDRH1_AA": CDRH1_AA, "CDRH3_AA": CDRH3_AA, "CDRL3_AA": CDRL3_AA}
            successful_updates += 1
        except Exception as e:
            print(f"Error updating sequence {id} (short: {id_short}): {e}")
            failed_updates += 1
    
    print(f"Dict_update results: {successful_updates} successful, {failed_updates} failed")
    return Seq_dict_update

def Check_quality(Seq_dict):
    L = [len(Seq_dict[i]["AA"]) for i in Seq_dict.keys()]
    print(f"Max:{max(L)}, Min:{min(L)}, Median: {np.median(L)}")
    return None

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run overlap region selection for antibody library construction")
    parser.add_argument('-n', '--negative', required=True, help='Path to negative control CSV file')
    parser.add_argument('-g', '--group-size', type=int, help='Group size from Step 4 (if not provided, will try to read from step4/group_size.txt)')
    parser.add_argument("-i", "--input", required=True, help="Path to input CSV file (same as used in Step 1)")
    
    args = parser.parse_args()
    
    # Get group size
    group_size = args.group_size
    if group_size is None:
        # Try to read from file saved by Step 4
        try:
            with open('ui_results/step4/group_size.txt', 'r') as f:
                group_size = int(f.read().strip())
        except (FileNotFoundError, ValueError):
            # Fallback to default
            group_size = 25
            print(f"Warning: Could not read group size from Step 4, using default: {group_size}")
    
    print(f"Using group_size: {group_size}")
    
    # Update the negative file path
    negative_file_path = args.negative
    
    # Main execution logic
    Group = [i for i in os.listdir("ui_results/step4") if i.startswith("Re_assembled")]
    TB = pd.read_csv(args.input, header=1)
    AA_L = 10
    os.system('mkdir -p ui_results/step5/primers')
    os.system('mkdir -p blastDB')
    os.system('mkdir -p ui_results/step5/blast')
    os.system('mkdir -p ui_results/step5/segs_id')

    for group in Group:
        print(f"Processing group: {group}")
        Fasta_DB = "ui_results/step4/" + group
        Seq_dict = read_F2D(Fasta_DB) # read the fasta
        Seq_dict = Dict_update(Seq_dict, TB) #Update Seq_dict by adding the seq from the table
        
        # Generate primers
        with open(f"ui_results/step5/primers/{group}", 'w') as F:
            for id in Seq_dict.keys():
                try:
                    i_CDRL3_AA = Seq_dict[id]["AA"].find(Seq_dict[id]["CDRL3_AA"]) #65
                except:
                    i_CDRL3_AA = 88
                i_CDRH1_AA = Seq_dict[id]["AA"].find(Seq_dict[id]["CDRH1_AA"]) #120
                i_CDRH3_AA = Seq_dict[id]["AA"].find(Seq_dict[id]["CDRH3_AA"]) #180
                for i in range(30):
                    Seq_fragment1 = Seq_dict[id]["Nr"][           :(i_CDRL3_AA + AA_L)*3-20+i]
                    Seq_fragment2 = Seq_dict[id]["Nr"][i_CDRL3_AA*3  :(i_CDRH1_AA + AA_L)*3-20+i]
                    Seq_fragment3 = Seq_dict[id]["Nr"][i_CDRH1_AA*3  :(i_CDRH3_AA + AA_L)*3-20+i]
                    F.write(f">{id}:Overlap1:{i}\n{Seq_fragment1[-30:]}\n")
                    F.write(f">{id}:Overlap2:{i}\n{Seq_fragment2[-30:]}\n")
                    F.write(f">{id}:Overlap3:{i}\n{Seq_fragment3[-30:]}\n")
        
        # Create BLAST database
        os.system(f"makeblastdb -in ui_results/step4/{group} -dbtype nucl -parse_seqids -out blastDB/{group}")
        
        # Initialize variables for BLAST loop
        TMP_tmp = pd.DataFrame()
        i = 0
        
        # Use configurable group_size instead of hardcoded 25
        target_sequences = group_size * 3
        print(f"Processing {group}: looking for {target_sequences} sequences (group_size={group_size} * 3)")
        
        while len(TMP_tmp) != target_sequences:
            i += 1
            # Run BLAST and process results
            os.system(f"blastn -query ui_results/step5/primers/{group} -db blastDB/{group} -out ui_results/step5/blast/{group} -outfmt '6 qseqid sseqid qstart qend sstart send' -num_threads 8 -max_hsps 2 -word_size {i}")
            os.system("cat ui_results/step5/blast/" + group + "| awk '{print $1}'| sort| uniq -c| awk '{print $1,$2}'| grep '^1 ' >  ui_results/step5/segs_id/" + group + '.list')
            
            # Read and process results inside the loop
            with open(f"ui_results/step5/segs_id/{group}.list", 'r') as F:
                Tmp = F.readlines()
            Tmp = [i.strip().split(" ")[1] for i in Tmp]
            TMP_tmp = pd.DataFrame({'uniq_id': Tmp})
            
            # Break if we've tried too many iterations to avoid infinite loop
            if i > 50:
                print(f"Warning: Maximum iterations reached for {group}, found {len(TMP_tmp)} sequences instead of {target_sequences}")
                break
        
        # Add metadata and save final results
        TMP_tmp['group'] = group
        TMP_tmp['word_size'] = i
        TMP_tmp.to_csv(f'ui_results/step5/segs_id/{group}.csv')
        print(f"Completed {group}: found {len(TMP_tmp)} sequences with word_size={i}")

if __name__ == "__main__":
    main()



