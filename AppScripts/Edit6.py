import streamlit as st
import pandas as pd
import os

# Streamlit App Title
st.title("🧬 Gene Expression Data Processor")

# Upload the two files
uploaded_file1 = st.file_uploader("📂 Upload IvyTumorInformation CSV", type="csv")
uploaded_file2 = st.file_uploader("📂 Upload Gene Expression CSV", type="csv")

if uploaded_file1 and uploaded_file2:
    # Load files into DataFrames
    file1 = pd.read_csv(uploaded_file1)
    file2 = pd.read_csv(uploaded_file2)

    # Show previews of uploaded files
    st.write("### 🔍 IvyTumorInformation Preview")
    st.dataframe(file1.head())

    st.write("### 🔍 Gene Expression Preview")
    st.dataframe(file2.head())

    # Process file2: Remove first column & transpose
    file2 = file2.iloc[:, 1:]
    file2_transposed = file2.transpose().reset_index()
    file2_transposed.columns = ['Gene'] + [f'Sample_{i}' for i in range(1, len(file2_transposed.columns))]

    # Ensure numeric values in gene expression columns
    for col in file2_transposed.columns[1:]:
        file2_transposed[col] = pd.to_numeric(file2_transposed[col], errors='coerce')

    # Select relevant columns from file1
    columns_to_select = ['donor_id', 'donor_name', 'sample_well', 'tumor_name', 'molecular_subtype', 'age_in_years']
    file1_selected = file1[columns_to_select]

    # Merge files based on 'sample_well' (file1) and 'Gene' (file2)
    merged_df = pd.merge(file1_selected, file2_transposed, left_on='sample_well', right_on='Gene', how='inner')

    # Handle potential data inconsistencies
    merged_df['Gene'] = merged_df['Gene'].apply(lambda x: '0' if str(x).count('.') == 2 else x)

    # Calculate average gene expression per donor
    gene_expression_data = merged_df.iloc[:, len(columns_to_select):]
    merged_df['AvgExpression'] = gene_expression_data.mean(axis=1)

    # Display the processed data
    st.write("### 📝 Processed Data Preview")
    st.dataframe(merged_df.head())

    # Save processed file
    processed_filename = "Processed_GeneExpression.csv"
    merged_df.to_csv(processed_filename, index=False, float_format="%.4f")

    # Download button for processed file
    st.download_button(
        label="📥 Download Processed Data",
        data=open(processed_filename, "rb"),
        file_name=processed_filename,
        mime="text/csv"
    )

    st.success("✅ Data processing complete! Download the file above.")

