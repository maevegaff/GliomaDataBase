import streamlit as st
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.formula.api import ols
import statsmodels.api as sm
from scipy.stats import shapiro, kruskal
import scikit_posthocs as sp
import itertools
import numpy as np

# Function to process & merge the files
def process_files(file1_path, file2_path):
    # Load files
    file1 = pd.read_csv(file1_path)
    file2 = pd.read_csv(file2_path)

    # Remove the first column in file2
    file2 = file2.iloc[:, 1:]

    # Transpose file2
    file2_transposed = file2.transpose().reset_index()

    # Rename columns
    file2_transposed.columns = ['Gene'] + [f'Sample_{i}' for i in range(1, len(file2_transposed.columns))]

    # Ensure numeric values
    for col in file2_transposed.columns[1:]:
        file2_transposed[col] = pd.to_numeric(file2_transposed[col], errors='coerce')

    # Select specific columns from file1
    columns_to_select = ['donor_id', 'donor_name', 'tumor_name', 'molecular_subtype']
    file1_selected = file1[columns_to_select]

    # Merge files
    merged_df = pd.concat([file1_selected, file2_transposed], axis=1)

    # Save processed data
    processed_file = "ProcessedData.csv"
    merged_df.to_csv(processed_file, index=False, float_format="%.4f")
    
    return processed_file, merged_df

# Function to generate and save boxplot with summary statistics
def generate_boxplot(data):
    plt.figure(figsize=(10, 6))
    ax = sns.boxplot(x="molecular_subtype", y="Gene", data=data, palette="Set2")
    
    # Perform Dunn's post hoc test
    dunn_result = sp.posthoc_dunn(data, val_col="Gene", group_col="molecular_subtype", p_adjust="bonferroni")

    # Get unique molecular subtypes
    subtypes = data["molecular_subtype"].unique()

    # Define y-position for significance markers
    y_max = data["Gene"].max()  
    y_offset = (y_max - data["Gene"].min()) * 0.05  
    y_pos = y_max + y_offset  

    # Iterate over all subtype pairs and annotate significance
    for (i, j) in itertools.combinations(range(len(subtypes)), 2):
        p_value = dunn_result.iloc[i, j]  
        if p_value < 0.05:  # Significant
            x1, x2 = i, j  
            ax.plot([x1, x1, x2, x2], [y_pos, y_pos + y_offset, y_pos + y_offset, y_pos], color="black", linewidth=1)
            ax.text((x1 + x2) / 2, y_pos + y_offset * 1.2, "*", ha="center", va="bottom", fontsize=14, color="red")
            y_pos += y_offset * 1.5  

    plt.xticks(rotation=45)
    plt.xlabel("Molecular Subtype")
    plt.ylabel("Gene Expression")
    plt.title("5HT3A Across Molecular Subtypes")

    # Save and display the plot
    plot_filename = "boxplot.png"
    plt.savefig(plot_filename, bbox_inches="tight")
    st.pyplot(plt)  # Display in Streamlit

    # Add a download button for the boxplot
    with open(plot_filename, "rb") as f:
        st.download_button(
            label="📥 Download Boxplot",
            data=f,
            file_name="Gene_Expression_Boxplot.png",
            mime="image/png"
        )
    
    return plot_filename

# Function to display summary statistics per molecular subtype
def display_summary_statistics(data):
    # Convert 'Gene' column to numeric (handling non-numeric values)
    data["Gene"] = pd.to_numeric(data["Gene"], errors="coerce")

    # Drop rows where 'Gene' is NaN
    data = data.dropna(subset=["Gene"])

    # Compute summary statistics per molecular subtype
    summary_stats = data.groupby("molecular_subtype")["Gene"].agg(
        count="count",
        mean="mean",
        std="std",
        min="min",
        q25=lambda x: x.quantile(0.25),
        median="median",
        q75=lambda x: x.quantile(0.75),
        max="max"
    )

    # Compute IQR (Interquartile Range)
    summary_stats["IQR"] = summary_stats["q75"] - summary_stats["q25"]

    # Rename columns for clarity
    summary_stats = summary_stats.rename(columns={"q25": "25%", "q75": "75%", "median": "Median", "std": "Std Dev"})

    st.write("### 📊 Molecular Subtype Summary Statistics")
    st.dataframe(summary_stats)

    # Download button for summary statistics
    summary_stats_csv = summary_stats.to_csv().encode('utf-8')
    st.download_button(
        label="📥 Download Summary Statistics",
        data=summary_stats_csv,
        file_name="MolecularSubtype_Summary_Statistics.csv",
        mime="text/csv"
    )

# Streamlit UI
st.title('🧬 Molecular Subtype Expression Analysis')

# Upload files
uploaded_file1 = st.file_uploader("📂 Upload IvyTumorInformation CSV", type="csv")
uploaded_file2 = st.file_uploader("📂 Upload Gene Expression CSV", type="csv")

if uploaded_file1 and uploaded_file2:
    # Save uploaded files
    file1_path = "uploaded_file1.csv"
    file2_path = "uploaded_file2.csv"
    
    with open(file1_path, "wb") as f:
        f.write(uploaded_file1.getbuffer())
    
    with open(file2_path, "wb") as f:
        f.write(uploaded_file2.getbuffer())

    # Process files
    processed_file, data = process_files(file1_path, file2_path)

    st.success("✅ Files merged successfully!")

    # Display preview
    st.write("### 🔍 Data Preview")
    st.dataframe(data.head())

    # Display summary statistics
    display_summary_statistics(data)

    # Generate and display boxplot
    st.write("### 📊 Gene Expression Boxplot")
    boxplot_file = generate_boxplot(data)

