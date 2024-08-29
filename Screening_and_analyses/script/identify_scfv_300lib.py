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
    df["Name"] = "" # Initialize the column outside the loop
    for index, row in df.iterrows():
        id += 1
        lev_dist = 100000
        for index_1, row in ref_df.iterrows():
            calculated_dist = Levenshtein.distance(str(df['muts'][index]), str(ref_df['nuc_seq'][index_1]))
            if calculated_dist < lev_dist:
                lev_dist = calculated_dist
                closest_abs = ref_df['Name'][index_1]
        if lev_dist < 1:
            df.loc[index, "Name"] = closest_abs
        else:
            df.loc[index, "Name"] = None
        #print(id)
    
    df = df.dropna(subset=["Name"])  # Drop rows with missing values in "Name" column
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

def get_score(df):
    df['Rep1_H1_stem_enrich'] = df['Rep1_H1_stem_postscreen_lib_freq'] / df['Rep1_prescreen_lib_freq'] 
    df['Rep2_H1_stem_enrich'] = df['Rep2_H1_stem_postscreen_lib_freq'] / df['Rep2_prescreen_lib_freq']
    df['Rep1_H3_stem_enrich'] = df['Rep1_H3_stem_postscreen_lib_freq'] / df['Rep1_prescreen_lib_freq'] 
    df['Rep2_H3_stem_enrich'] = df['Rep2_H3_stem_postscreen_lib_freq'] / df['Rep2_prescreen_lib_freq']
    return df

def get_avg_score(df):
    df['H1_stem_avg_enrich'] = (df['Rep1_H1_stem_enrich'] + df['Rep2_H1_stem_enrich'])/2
    df['H3_stem_avg_enrich'] = (df['Rep1_H3_stem_enrich'] + df['Rep2_H3_stem_enrich'])/2
    return df

def process_count_data(infile, ref_file, info_file, outfile_raw, outfile_info, num_processes, chunk_size):
    ref_df = read_tsv(ref_file)
    df = read_tsv(infile)
    info_df = pd.read_csv(info_file)
    print(info_df.columns)

    chunks = chunk_dataframe(df, chunk_size)
    pool = Pool(processes=num_processes)

    print("start mulitparallel processing")
    results = pool.starmap(identify_scFv, [(chunk, ref_df) for chunk in chunks])
    pool.close()
    pool.join()

    print("finish mulitparallel processing")

    combined_result_df = pd.concat(results)
    combined_result_df = get_freq(combined_result_df)
    combined_result_df = get_score(combined_result_df)
    combined_result_df = get_avg_score(combined_result_df)
    col_raw = ['Name', 'muts','Rep1_prescreen_lib_count', 'Rep2_prescreen_lib_count', 'Rep1_H1_stem_postscreen_lib_count', 
            'Rep2_H1_stem_postscreen_lib_count', 'Rep1_H3_stem_postscreen_lib_count', 'Rep2_H3_stem_postscreen_lib_count',
            'Rep1_prescreen_lib_freq', 'Rep2_prescreen_lib_freq', 'Rep1_H1_stem_postscreen_lib_freq', 
            'Rep2_H1_stem_postscreen_lib_freq', 'Rep1_H3_stem_postscreen_lib_freq', 'Rep2_H3_stem_postscreen_lib_freq',
            'Rep1_H1_stem_enrich','Rep2_H1_stem_enrich','Rep1_H3_stem_enrich','Rep2_H3_stem_enrich',
            'H1_stem_avg_enrich', 'H3_stem_avg_enrich']
    combined_result_df = combined_result_df[col_raw]
    combined_result_df.to_csv(outfile_raw, sep='\t', index=False)   

    filtered_info_df = info_df[['Name','Heavy_V_gene','Heavy_J_gene','Heavy_D_gene','Light_V_gene','Light_J_gene', 'Binds to','Reference']]

    
    combined_result_df = combined_result_df.merge(filtered_info_df , on= "Name", how = "inner")

    col_info = ['Name', 'muts','Rep1_prescreen_lib_count', 'Rep2_prescreen_lib_count', 'Rep1_H1_stem_postscreen_lib_count', 
            'Rep2_H1_stem_postscreen_lib_count', 'Rep1_H3_stem_postscreen_lib_count', 'Rep2_H3_stem_postscreen_lib_count',
            'Rep1_prescreen_lib_freq', 'Rep2_prescreen_lib_freq', 'Rep1_H1_stem_postscreen_lib_freq', 
            'Rep2_H1_stem_postscreen_lib_freq', 'Rep1_H3_stem_postscreen_lib_freq', 'Rep2_H3_stem_postscreen_lib_freq',
            'Rep1_H1_stem_enrich','Rep2_H1_stem_enrich','Rep1_H3_stem_enrich','Rep2_H3_stem_enrich',
            'H1_stem_avg_enrich', 'H3_stem_avg_enrich','Heavy_V_gene','Heavy_J_gene','Heavy_D_gene','Light_V_gene','Light_J_gene', 'Binds to','Reference']
    combined_result_df = combined_result_df[col_info]
    combined_result_df.to_csv(outfile_info, sep='\t', index=False)    

    print("finish score calculation")

    return combined_result_df

def main():
    start_time = time.time()
    ref_file = 'ref_files/300lib.tsv'
    info_file = 'ref_files/300lib_Abs.csv'
    outfile_raw = "result/PacBio/oPool_screen_counts_freq_and_enrichment.tsv"
    outfile_info = "result/PacBio/oPool_screen_enrichment.tsv"
    infile = 'result/PacBio/mut_nuc_count.tsv'

    total_time = time.time() - start_time

    process_count_data(infile, ref_file, info_file, outfile_raw, outfile_info, 10, 20000)

    total_time = time.time() - start_time
    print(f"Total processing time: {total_time:.2f} seconds")

if __name__ == "__main__":
    main()
