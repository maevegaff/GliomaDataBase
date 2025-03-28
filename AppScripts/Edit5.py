import pandas as pd
import os
import requests


# Import the two files
#file1 consists of patient and tumor details
#file2 is gene expression file

#This file path should be according to the local path download 
file_path1 = r"C:\Users\maeve\Downloads\IvyTumorInfomation.csv"
file1 = pd.read_csv(file_path1)

print(file1.head())

#This file path should be according to the local path download 
file_path2 = r"C:\Users\maeve\Downloads\Expression5HT3A.csv"
file2 = pd.read_csv(file_path2)

print(file2.head())

# Delete the first column in file2
file2 = file2.iloc[:, 1:]

# Transpose the rows of file2 into columns
file2_transposed = file2.transpose()

# Reset the index to make the transposed data easier to work with
file2_transposed.reset_index(inplace=True)

# Rename the columns to reflect the transposition
file2_transposed.columns = ['Gene'] + [f'Sample_{i}' for i in range(1, len(file2_transposed.columns))]

# Ensure gene expression values are treated as floats
for col in file2_transposed.columns[1:]:
    file2_transposed[col] = file2_transposed[col].astype(float)


# Display the transposed dataframe
print("\nTransposed file2:")
print(file2_transposed)

# Create new dataframe with selected columns from file1
columns_to_select_file1 = ['donor_id', 'donor_name', 'donor_color', 'sample_well', 'sample_polygon', 'sample_mri_0', 'sample_mri_1', 'sample_mri_2', 'structure_id', 'structure_name', 'structure_abbreviation', 'structure_color', 'top_level_structure_id', 'top_level_structure_abbreviation', 'top_level_structure_color', 'tumor_name', 'molecular_subtype', 'extent_of_resection', 
    'surgery', 'mgmt_methylation', 'survival_days', 'egfr_amplification', 
    'initial_kps', 'age_in_years'] # These are the columns I want to select from file1

file1_selected = file1[columns_to_select_file1]  # This correctly selects multiple columns

# Merge the selected columns from file1 with the transposed file2
# Assuming 'sample_well' in file1_selected corresponds to 'Gene' in file2_transposed
new_df = pd.merge(file1_selected, file2_transposed, left_on='sample_well', right_on='Gene', how='inner')
new_df['Gene'] = new_df['Gene'].apply(lambda x: '0' if str(x).count('.') == 2 else x)
new_df['Gene'] = new_df['Gene'].apply(lambda x: '0' if x.count('.') == 2 else x)

# Calculate the mean gene expression for each donor_id
# Extract gene expression data
gene_expression_data = new_df.iloc[:, 24:]

# Handle NaN values in donor_id
new_df['donor_id'].fillna('Unknown', inplace=True)

# Calculate the overall mean expression
new_df['AvgExpression'] = mean_expression_per_donor.mean(axis=0)

# Calculate the overall mean expression
new_df['AvgExpression'] = mean_expression_per_donor.mean(axis=1)

# Display the dataframe with the new column
print("\nDataFrame with AvgExpression column:")
print(new_df.head())


# Save the new file
new_df.to_csv('new_file.csv', index=False, float_format='%.4f')








