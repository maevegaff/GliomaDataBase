import pandas as pd
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Merge tumor and gene expression data.")
parser.add_argument("--file1", required=True, help="Path to the tumor information CSV file")
parser.add_argument("--file2", required=True, help="Path to the gene expression CSV file")
parser.add_argument("--output", required=True, help="Path to the output CSV file")

args = parser.parse_args()

# Load the input files
file1 = pd.read_csv(args.file1)
file2 = pd.read_csv(args.file2)

# Delete the first column in file2
file2 = file2.iloc[:, 1:]

# Transpose file2
file2_transposed = file2.transpose()
file2_transposed.reset_index(inplace=True)
file2_transposed.columns = ['Gene'] + [f'Sample_{i}' for i in range(1, len(file2_transposed.columns))]

# Ensure gene expression values are treated as floats
for col in file2_transposed.columns[1:]:
    file2_transposed[col] = file2_transposed[col].astype(float)

# Select specific columns from file1
columns_to_select = ['donor_id', 'donor_name', 'donor_color', 'sample_well', 'sample_polygon', 'sample_mri_0', 
                     'sample_mri_1', 'sample_mri_2', 'structure_id', 'structure_name', 'structure_abbreviation', 
                     'structure_color', 'top_level_structure_id', 'top_level_structure_abbreviation', 
                     'top_level_structure_color', 'tumor_name', 'molecular_subtype', 'extent_of_resection', 
                     'surgery', 'mgmt_methylation', 'survival_days', 'egfr_amplification', 'initial_kps', 'age_in_years']

file1_selected = file1[columns_to_select]

# Merge the data
new_df = pd.concat([file1_selected, file2_transposed], axis=1)

# Adjust 'Gene' column formatting
new_df['Gene'] = new_df['Gene'].apply(lambda x: '0' if x.count('.') == 2 else x)
# Save output
new_df.to_csv(args.output, index=False, float_format='%.4f')

# Save the new data frame under the name GalaxyExpression
GalaxyExpression = new_df
