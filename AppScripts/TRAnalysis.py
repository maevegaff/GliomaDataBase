import streamlit as st
import pandas as pd
import os
import glob
from statsmodels.formula.api import ols
import statsmodels.api as sm
from scipy.stats import shapiro, kruskal
import scikit_posthocs as sp

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
    

    # TumorRegionAnalysis

    # Load the data
    data = pd.read_csv(new_file_path)

    # Print the data
    st.write(data.head())

    # Summary statistics
    summary_stats = data.describe()
    st.write(summary_stats)

    # Save the summary statistics to a CSV file
    summary_stats_file = 'SumStat.csv'
    summary_stats.to_csv(summary_stats_file)
    st.download_button(
        label="Download Summary Statistics",
        data=open(summary_stats_file, "rb").read(),
        file_name=summary_stats_file,
        mime="text/csv"
    )

    # Define the model
    model = ols('Gene ~ C(structure_color)', data=data).fit()

    # Perform ANOVA
    anova_table = sm.stats.anova_lm(model, typ=2)
    st.write(anova_table)

    # Save the ANOVA results to a CSV file
    anova_results_file = 'ANOVAResults.csv'
    anova_table.to_csv(anova_results_file)
    st.download_button(
        label="Download ANOVA Results",
        data=open(anova_results_file, "rb").read(),
        file_name=anova_results_file,
        mime="text/csv"
    )

    # Shapiro-Wilk Test
    stat, p = shapiro(data['Gene'])  # Replace 'Gene' with actual column name if different
    st.write(f"Shapiro-Wilk Test: p-value = {p}")

    # Run Kruskal-Wallis Test
    stat, p = kruskal(*[data[data['structure_color'] == region]['Gene'] for region in data['structure_color'].unique()])
    st.write(f"Kruskal-Wallis Test: p = {p}")

    # If significant, run Dunn’s test
    dunn = sp.posthoc_dunn(data, val_col='Gene', group_col='structure_color', p_adjust='bonferroni')
    st.write(dunn)

    # Save the results to a CSV file
    dunn_results_file = 'DunnResults.csv'
    dunn.to_csv(dunn_results_file)
    st.download_button(
        label="Download Dunn's Test Results",
        data=open(dunn_results_file, "rb").read(),
        file_name=dunn_results_file,
        mime="text/csv"
    )

import streamlit as st
import pandas as pd
import os
from statsmodels.formula.api import ols
import statsmodels.api as sm
from scipy.stats import shapiro, kruskal
import scikit_posthocs as sp

def process_file(file_path):
    # Load the file
    data = pd.read_csv(file_path)

    # Define the model
    model = ols('Gene ~ C(structure_color)', data=data).fit()

    # Perform ANOVA
    anova_table = sm.stats.anova_lm(model, typ=2)
    st.write(anova_table)

    # Save the ANOVA results to a CSV file
    anova_results_file = 'ANOVAResults.csv'
    anova_table.to_csv(anova_results_file)
    st.download_button(
        label="Download ANOVA Results",
        data=open(anova_results_file, "rb").read(),
        file_name=anova_results_file,
        mime="text/csv"
    )

    # Shapiro-Wilk Test
    stat, p = shapiro(data['Gene'])  # Replace 'Gene' with actual column name if different
    st.write(f"Shapiro-Wilk Test: p-value = {p}")

    # Run Kruskal-Wallis Test
    stat, p = kruskal(*[data[data['structure_color'] == region]['Gene'] for region in data['structure_color'].unique()])
    st.write(f"Kruskal-Wallis Test: p = {p}")

    # If significant, run Dunn’s test
    dunn = sp.posthoc_dunn(data, val_col='Gene', group_col='structure_color', p_adjust='bonferroni')
    st.write(dunn)

    # Save the results to a CSV file
    dunn_results_file = 'DunnResults.csv'
    dunn.to_csv(dunn_results_file)
    
    # Display the link to download the Dunn's test results
    st.success(f"Dunn's test results processed successfully! [Download DunnResults.csv](./{dunn_results_file})")
    st.download_button(
        label="Download Dunn's Test Results",
        data=open(dunn_results_file, "rb").read(),
        file_name="DunnResults.csv",
        mime="text/csv"
    )

st.title('Tumour Region Analysis')

# Upload the new_file.csv
uploaded_file = st.file_uploader("Choose the new_file.csv", type="csv")

if uploaded_file is not None:
    # Save the uploaded file to disk
    file_path = os.path.join("new_file.csv")
    
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    
    # Process the file
    process_file(file_path)
