import pandas as pd


input_file = 'oPool_result/202412_mut_aa_count.tsv'  
df = pd.read_csv(input_file, sep='\t')

file_names_df = pd.read_csv('ref_files/202412_sample_name.tsv', sep='\t')
file_names_dict = file_names_df.set_index('sample_ID')['sample_name'].to_dict()

subset_columns = {
    'Full_HA_CR9114_competition': ['muts'] + [str(i) for i in range(1, 15)],  
    'Full_HA': ['muts'] + [str(i) for i in range(15, 31)], 
    'Assembly_25': ['muts'] + ['31', '37'], 
    'Assembly_50': ['muts'] + ['32', '38'], 
    'Assembly_75': ['muts'] + ['33', '39'], 
    'Assembly_100': ['muts'] + ['34', '40'], 
    'Assembly_125': ['muts'] + ['35', '41'], 
    'Assembly_150': ['muts'] + ['36', '42'],
    'Assembly_200': ['muts'] + ['45', '46']
}

print(subset_columns)


for subset_name, columns in subset_columns.items():
    subset_df = df[columns] 
    numerical_columns = subset_df.select_dtypes(include=[int, float]).columns
    df_filtered = subset_df[(subset_df[numerical_columns] != 0).any(axis=1)]
    for column_name in subset_df.columns:
        if column_name != "muts":
            new_column_name = file_names_dict.get(int(column_name), f"Column_{column_name}") + '_count'
            df_filtered = df_filtered.rename(columns={column_name: new_column_name})

    df_filtered.to_csv(f'oPool_result/nuc_count_files/{subset_name}.tsv', sep='\t', index=False) 
