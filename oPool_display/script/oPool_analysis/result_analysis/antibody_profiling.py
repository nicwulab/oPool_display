import pandas as pd
import os

# List of file paths
file_paths = [
    'result/filtered_antigen_Lee_FluB_custom_cutoff.tsv',
    'result/filtered_antigen_Phu_FluB_custom_cutoff.tsv',
    'result/filtered_antigen_SH_H7_custom_cutoff.tsv',
    'result/filtered_antigen_QH_H5_custom_cutoff.tsv',
    'result/filtered_antigen_SP16_H3_custom_cutoff.tsv',
    'result/filtered_antigen_MI15_H1_custom_cutoff.tsv',
    'result/filtered_antigen_SI06_H1_custom_cutoff.tsv',
    'result/filtered_antigen_H3_stem_custom_cutoff.tsv',
    'result/filtered_antigen_H1_stem_custom_cutoff.tsv'
]

# Mapping file names to antigen types
file_to_antigen = {
    'filtered_antigen_Lee_FluB_custom_cutoff.tsv': 'B/Lee40',
    'filtered_antigen_Phu_FluB_custom_cutoff.tsv': 'B/Phu13',
    'filtered_antigen_SH_H7_custom_cutoff.tsv': 'H7/SH13',
    'filtered_antigen_QH_H5_custom_cutoff.tsv': 'H5/QH05',
    'filtered_antigen_SP16_H3_custom_cutoff.tsv': 'H3/SP16',
    'filtered_antigen_MI15_H1_custom_cutoff.tsv': 'H1/MI15',
    'filtered_antigen_SI06_H1_custom_cutoff.tsv': 'H1/SI06',
    'filtered_antigen_H3_stem_custom_cutoff.tsv': 'H3 stem',
    'filtered_antigen_H1_stem_custom_cutoff.tsv': 'H1 stem',
}

# Define dynamic competition index cutoff values for each antigen
comp_index_cutoffs = {
    'B/Lee40': 1.4, 
    'B/Phu13': 1.65, 
    'H7/SH13': 2,
    'H5/QH05': 2.5,
    'H1/MI15': 2.5,
    'H1/SI06': 2
}

# Load antibody metadata from 300lib_Abs.tsv
metadata_path = 'ref_files/300lib_Abs.tsv'
metadata_df = pd.read_csv(metadata_path, sep='\t')

# Select relevant columns
selected_columns = ['VH_AA', 'VL_AA', 'PMID']
selected_columns += [col for col in metadata_df.columns if '_gene' in col]
rename_mapping = {
    'Specificity': 'Prior knowledge: specificity',
    'Binds to': 'Prior knowledge: binds to'
}
metadata_df.rename(columns=rename_mapping, inplace=True)
selected_columns += list(rename_mapping.values())

metadata_df = metadata_df[['Name'] + selected_columns]
metadata_df.set_index('Name', inplace=True)

# Dictionary to store antibody information
antibody_info = {}

# Reading each file and collecting antigen specificity
for file in file_paths:
    filename = os.path.basename(file).strip()
    if filename not in file_to_antigen:
        print(f"Warning: Filename '{filename}' not found in file_to_antigen mapping.")
        continue

    antigen = file_to_antigen[filename]
    comp_index_cutoff = comp_index_cutoffs.get(antigen, 2.5)  # Default to 2.5 if not found

    df = pd.read_csv(file, sep='\t')
    if 'closest_abs' in df.columns:
        for _, row in df.iterrows():
            ab = row['closest_abs']
            if pd.notna(ab):
                if ab not in antibody_info:
                    antibody_info[ab] = {
                        'specificity': set(), 
                        'number_of_antigen': 0,
                        'compete_with_CR9114_on': set(),
                        'do_not_compete_with_CR9114_on': set()
                    }
                antibody_info[ab]['specificity'].add(antigen)
                antibody_info[ab]['number_of_antigen'] += 1

                # Check for competition index with dynamic cutoff
                if 'competition_index' in df.columns:
                    comp_index = row['competition_index']
                    if pd.notna(comp_index):
                        if comp_index > comp_index_cutoff:
                            antibody_info[ab]['compete_with_CR9114_on'].add(antigen)
                        else:
                            antibody_info[ab]['do_not_compete_with_CR9114_on'].add(antigen)

# Preparing the final DataFrame
antibody_list = []
for ab, info in antibody_info.items():
    row_data = {
        'Name': ab,
        'specificity': ', '.join(sorted(info['specificity'])) if info['specificity'] else 'N/A',
        'number_of_antigen': info['number_of_antigen'],
        'compete_with_CR9114_on': ', '.join(sorted(info['compete_with_CR9114_on'])) if info['compete_with_CR9114_on'] else 'N/A',
        'do_not_compete_with_CR9114_on': ', '.join(sorted(info['do_not_compete_with_CR9114_on'])) if info['do_not_compete_with_CR9114_on'] else 'N/A'
    }
    # Add metadata from 300lib_Abs.tsv
    if ab in metadata_df.index:
        row_data.update(metadata_df.loc[ab].to_dict())
    antibody_list.append(row_data)

# Convert to DataFrame
antibody_df = pd.DataFrame(antibody_list)

# Save to TSV file
output_path = 'result/unique_antibodies_info.tsv'
antibody_df.to_csv(output_path, sep='\t', index=False)

print(f"The output file has been saved as {output_path}")