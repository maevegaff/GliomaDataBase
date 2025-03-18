import streamlit as st
import pandas as pd
import os
import glob

def process_files(file1_path, file2_path):
    # Load the files
    file1 = pd.read_csv(file1_path)
    file2 = pd.read_csv(file2_path)

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

    # Create new dataframe with selected columns from file1
    columns_to_select_file1 = ['donor_id', 'donor_name', 'donor_color', 'sample_well', 'sample_polygon', 'sample_mri_0', 'sample_mri_1', 'sample_mri_2', 'structure_id', 'structure_name', 'structure_abbreviation', 'structure_color', 'top_level_structure_id', 'top_level_structure_abbreviation', 'top_level_structure_color', 'tumor_name', 'molecular_subtype', 'extent_of_resection', 
        'surgery', 'mgmt_methylation', 'survival_days', 'egfr_amplification', 
        'initial_kps', 'age_in_years'] # These are the columns I want to select from file1

    file1_selected = file1[columns_to_select_file1]  # This correctly selects multiple columns

    # Merge the selected columns from file1 with the transposed file2
    new_df = pd.concat([file1_selected, file2_transposed], axis=1)

    # Save the new file
    new_file_path = 'new_file.csv'
    new_df.to_csv(new_file_path, index=False, float_format='%.4f')

    return new_file_path

def run_workflow_scripts():
    script_dir = '/c:/Users/maeve/.vscode/GliomaDataBase/WorkflowScripts/'
    script_files = glob.glob(os.path.join(script_dir, '*.py'))

    for script_file in script_files:
        with open(script_file) as f:
            code = compile(f.read(), script_file, 'exec')
            exec(code, globals())

st.title('Gene Expression Analysis')

# Upload the gene expression file
uploaded_file1 = st.file_uploader("Choose the IvyTumorInformation CSV file", type="csv")
uploaded_file2 = st.file_uploader("Choose the Gene Expression CSV file", type="csv")


if uploaded_file1 is not None and uploaded_file2 is not None:
    # Save the uploaded files to disk
    file1_path = os.path.join("uploaded_file1.csv")
    file2_path = os.path.join("uploaded_file2.csv")
    
    with open(file1_path, "wb") as f:
        f.write(uploaded_file1.getbuffer())
    
    with open(file2_path, "wb") as f:
        f.write(uploaded_file2.getbuffer())
    
    # Process the files
    new_file_path = process_files(file1_path, file2_path)
    
    # Display the link to download the new file
    st.success(f"File processed successfully! [Download new_file.csv](./{new_file_path})")