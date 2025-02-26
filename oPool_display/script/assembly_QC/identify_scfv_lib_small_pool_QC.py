import pandas as pd
from Bio import SeqIO
import csv
from collections import Counter
from functools import reduce
import operator
import Levenshtein
import time
from multiprocessing import Pool

def read_tsv(file_path):
    df = pd.read_csv(file_path, delimiter='\t')
    return df



def calculate_levenshtein_distance(string1, string2):
    distance = Levenshtein.distance(string1, string2)
    return distance

def create_ref_file(list_of_tsv_files):
    ref_df = pd.DataFrame(columns=['name', 'nuc_seq'])
    for tsv_file in list_of_tsv_files:
        df = read_tsv(tsv_file)
        ref_df = pd.concat([ref_df, df])
    ref_df = ref_df.drop_duplicates()
    ref_df.reset_index(drop=True, inplace=True)
    ref_df.to_csv('ref_seq.tsv', sep='\t', index=False)
    print(ref_df)
    return ref_df

def identify_scFv(df, ref_df):
   # print("Start identify scFv..")
    id = 0 
    df["closest_abs"] = "" # Initialize the column outside the loop
    df["lev_dist"] = ""
    for index, row in df.iterrows():
        id += 1
        lev_dist = 100000
        for index_1, row in ref_df.iterrows():
            calculated_dist = Levenshtein.distance(str(df['muts'][index]), str(ref_df['nuc_seq'][index_1]))
            if calculated_dist < lev_dist:
                lev_dist = calculated_dist
                closest_abs = ref_df['name'][index_1]
        if lev_dist < 1:
            df.loc[index, "closest_abs"] = closest_abs
            df.loc[index, "lev_dist"] = lev_dist
        else:
            df.loc[index, "closest_abs"] = None
        #print(id)
    
    df = df.dropna(subset=["closest_abs"])  # Drop rows with missing values in "closest_abs" column
   # print("Done")
    return df

def chunk_dataframe(df, chunk_size):
    chunks = [df[i:i+chunk_size] for i in range(0, len(df), chunk_size)]
    return chunks

def count_to_freq_col(df, colname):
    df[colname] = pd.to_numeric(df[colname], errors='coerce')
    new_col_name = colname[:-6] + '_freq'
    print('calculate freq for: ' + colname[:-6])
    df[new_col_name] = (df[colname] + 1) / (df[colname].sum() + len(df))
    return df

def apply_count_to_freq_parallel(df, columns):
    with Pool() as pool:
        results = pool.starmap(count_to_freq_col, [(df, col) for col in columns])
    return results

def count_to_freq(df, colname):
    df[colname[:-6]+'_freq'] = (df[colname]+1)/(df[colname].sum()+len(df))
    return df

def get_freq(df):
    colnames = [colname for colname in df.columns]
    for col in colnames:
        if 'count' in col:
            count_to_freq(df, col)
    return df


def process_count_data(infile, ref_file, outfile, num_processes, chunk_size):
    ref_df = read_tsv(ref_file)
    df = read_tsv(infile)

    chunks = chunk_dataframe(df, chunk_size)
    pool = Pool(processes=num_processes)

    print("start mulitparallel processing")
    results = pool.starmap(identify_scFv, [(chunk, ref_df) for chunk in chunks])
    pool.close()
    pool.join()

    print("finish mulitparallel processing")

    combined_result_df = pd.concat(results)
    combined_result_df = get_freq(combined_result_df)
    combined_result_df.to_csv(outfile, sep='\t', index=False)

    print("finish score calculation")

    return combined_result_df

def filter_ref_file(ref_file, name_file):
    ref_df = read_tsv(ref_file)
    name_df = pd.read_csv(name_file, header=None, names=["name"])
    filtered_ref_df = ref_df[ref_df['name'].isin(name_df['name'])]
    print(filtered_ref_df)
    return filtered_ref_df

def main():
    start_time = time.time()

    # Single reference file with all sequences
    ref_file = 'ref_files/lib_ref.csv'

    # List of name-only files
    name_files = [
        "ref_files/small_pool_lists/25.tsv",
        "ref_files/small_pool_lists/50.tsv",
        "ref_files/small_pool_lists/75.tsv",
        "ref_files/small_pool_lists/100.tsv",
        "ref_files/small_pool_lists/125.tsv",
        "ref_files/small_pool_lists/150.tsv",
        "ref_files/small_pool_lists/200.tsv"
    ]

    # List of input data files
    infile_list = [
        "oPool_result/nuc_count_files/Assembly_25.tsv",
        "oPool_result/nuc_count_files/Assembly_50.tsv",
        "oPool_result/nuc_count_files/Assembly_75.tsv",
        "oPool_result/nuc_count_files/Assembly_100.tsv",
        "oPool_result/nuc_count_files/Assembly_125.tsv",
        "oPool_result/nuc_count_files/Assembly_150.tsv",
        "oPool_result/nuc_count_files/Assembly_200.tsv"
    ]

    # Corresponding output files
    outfile_list = [
        "oPool_result/assembly_QC/freq/25_assembly_freq.tsv",
        "oPool_result/assembly_QC/freq/50_assembly_freq.tsv",
        "oPool_result/assembly_QC/freq/75_assembly_freq.tsv",
        "oPool_result/assembly_QC/freq/100_assembly_freq.tsv",
        "oPool_result/assembly_QC/freq/125_assembly_freq.tsv",
        "oPool_result/assembly_QC/freq/150_assembly_freq.tsv",
        "oPool_result/assembly_QC/freq/200_assembly_freq.tsv"
    ]

    for name_file, infile, outfile in zip(name_files, infile_list, outfile_list):
        # Filter reference file based on the current name file
        filtered_ref_df = filter_ref_file(ref_file, name_file)

        # Process the input file using the filtered reference DataFrame
        process_count_data(infile, filtered_ref_df, outfile, 60, 1000)

    total_time = time.time() - start_time
    print(f"Total processing time: {total_time:.2f} seconds")

if __name__ == "__main__":
    main()