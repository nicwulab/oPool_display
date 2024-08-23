#!/usr/bin/python
import operator
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
import pandas as pd
import multiprocessing
import time
import os
import re

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

def translation(seq):
    dnamap = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
              "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
              "TAT": "Y", "TAC": "Y", "TAA": "_", "TAG": "_",
              "TGT": "C", "TGC": "C", "TGA": "_", "TGG": "W",
              "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
              "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
              "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
              "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
              "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
              "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
              "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
              "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
              "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
              "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
              "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
              "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", }
    pep = []
    i = 0
    while i < len(seq):
        codon = seq[i:i + 3]
        aa = dnamap[codon]
        pep.append(aa)
        i = i + 3
    pep = ''.join(pep)
    return pep


def extract_sequences_between_motifs(sequence, motif1, motif2):
    # Convert input sequences to strings
    seq_str = str(sequence)
    motif1_str = str(motif1)
    motif2_str = str(motif2)

    # Find matches for the forward motifs
    match = None
    for m1 in re.finditer(motif1_str, seq_str):
        for m2 in re.finditer(motif2_str, seq_str):
            if m1.end() < m2.start():
                match = seq_str[m1.end():m2.start()].upper()
                break
        if match:
            break
    
    # If no matches found, try reverse complement motifs
    if match is None:
        seq = Seq(sequence)
        rev_compl_seq = seq.reverse_complement()
        rev_compl_seq_str = str(rev_compl_seq)
        for m1 in re.finditer(str(Seq(motif1).reverse_complement()), rev_compl_seq_str):
            for m2 in re.finditer(str(Seq(motif2).reverse_complement()), rev_compl_seq_str):
                if m1.end() < m2.start():
                    match = rev_compl_seq_str[m1.end():m2.start()].upper()
                    break
            if match:
                break
    
    if match is not None:
        return match


def extract_sequences_between_motifs_from_fastq(fastq_file):
    motif1 = 'GGAGTATCCACCATG'
    motif2 = 'TCCGGAGGATCCGAT'
    matches = []
    with open(fastq_file, "r") as handle:
        record_count = 0
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            matches.append(extract_sequences_between_motifs(seq, motif1, motif2))
    return matches

def process_seq(fastq_file):
    forward_motif = 'GGAGTATCCACCATG'
    variants = []
    with open(fastq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            if seq[:15] == forward_motif:
                variants.append(seq[15:-15])
            else:
                reverse_complement_seq = Seq(seq).reverse_complement()
                variants.append(str(reverse_complement_seq)[15:-15])
    return variants

def ProcessMultilib(Rfile):
    print("Reading %s" % Rfile)
    variants = []
    variants =  process_seq(Rfile)
    return Counter(variants)


def process_file(fastq_file, file_names_dict):
    #reffile = 'fasta/ref_seq.fasta' 
    #refseq  = next(SeqIO.parse(reffile,"fasta")).seq
    count_dictionary = ProcessMultilib(fastq_file)

    split_filename = fastq_file.split("_")
    split_further = split_filename[1].split("/")
    sample_name = split_further[1]
    name_info = file_names_dict[int(sample_name)]
    #print(name_info)
    return sample_name, count_dictionary


def main():
    start_time = time.time()

    outfile = 'result/PacBio/mut_nuc_count.tsv'
    fastq_list = glob.glob('fastq_filtered/*.fastq')
    file_names_df = pd.read_csv('ref_files/sample_name.tsv', sep='\t')
    file_names_dict = file_names_df.set_index('sample_ID')['sample_name'].to_dict()

    count_df = pd.DataFrame()
    count = 0
    num_processes = 8

    count_df = None
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.starmap(process_file, [(fastq_file, file_names_dict) for fastq_file in fastq_list])

    for name_info, count_dictionary in results:
        if count_df is None:
            count_df = pd.DataFrame.from_records(list(dict(count_dictionary).items()), columns=['muts', name_info])
            count_df = count_df[count_df['muts'].notna()] 
            count_df['muts'] = count_df['muts'].apply(lambda tup: ''.join(map(str, tup)))
        else:
            sub_df = pd.DataFrame.from_records(list(dict(count_dictionary).items()), columns=['muts', name_info])
            sub_df = sub_df[sub_df['muts'].notna()] 
            sub_df['muts'] = sub_df['muts'].apply(lambda tup: ''.join(map(str, tup)))
            count_df = count_df.merge(sub_df, on='muts', how='outer')
    
    for column_name in count_df.columns:
        if column_name != "muts":
            new_column_name = file_names_dict.get(int(column_name), f"Column_{column_name}") + '_count'
            count_df = count_df.rename(columns={column_name: new_column_name})

    count_df = count_df.fillna(0)

    cols = ['muts', 'Rep1_prescreen_lib_count', 'Rep2_prescreen_lib_count', 'Rep1_H1_stem_postscreen_lib_count', 
            'Rep1_H3_stem_postscreen_lib_count', 'Rep2_H1_stem_postscreen_lib_count', 'Rep2_H3_stem_postscreen_lib_count']
  
    count_df = count_df[cols]
    print(count_df)
    count_df.to_csv(outfile, sep="\t", index = False)

    total_time = time.time() - start_time
    print(f"Total processing time: {total_time:.2f} seconds")

if __name__ == "__main__":
  main()
