import streamlit as st
import pandas as pd
import os
import glob
from statsmodels.formula.api import ols
import statsmodels.api as sm
from scipy.stats import shapiro, kruskal
import scikit_posthocs as sp

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
    columns_to_select = [
        'donor_id', 'donor_name', 'structure_color', 'tumor_name', 'molecular_subtype'
    ]  # Add more if needed
    file1_selected = file1[columns_to_select]

    # Merge files
    merged_df = pd.concat([file1_selected, file2_transposed], axis=1)

    # Save processed data
    processed_file = "ProcessedData.csv"
    merged_df.to_csv(processed_file, index=False, float_format="%.4f")
    
    return processed_file, merged_df

# Streamlit UI
st.title('🧬 Tumor Region Expression Analysis')

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

    # Download button for processed data
    st.success("✅ Files merged successfully!")
    st.download_button(
        label="📥 Download Processed Data",
        data=open(processed_file, "rb").read(),
        file_name="ProcessedData.csv",
        mime="text/csv"
    )

    # Display preview
    st.write("### 🔍 Data Preview")
    st.dataframe(data.head())

    # Ensure 'Gene' and 'structure_color' exist
    required_columns = ["structure_color", "Gene"]
    if not all(col in data.columns for col in required_columns):
        st.error(f"❌ CSV must contain: {required_columns}")
        st.stop()

    # Convert Gene column to numeric
    try:
        data["Gene"] = pd.to_numeric(data["Gene"], errors="coerce")
    except ValueError:
        st.error("❌ 'Gene' column must contain numeric values.")
        st.stop()

    # Summary statistics
    summary_stats = data.describe()
    st.write("### 📊 Summary Statistics")
    st.dataframe(summary_stats)

    # Save and download summary statistics
    summary_stats_file = "SummaryStatistics.csv"
    summary_stats.to_csv(summary_stats_file)
    st.download_button(
        label="📥 Download Summary Statistics",
        data=open(summary_stats_file, "rb").read(),
        file_name="SummaryStatistics.csv",
        mime="text/csv"
    )

    # Perform ANOVA
    model = ols("Gene ~ C(structure_color)", data=data).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)

    st.write("### 📊 ANOVA Results")
    st.dataframe(anova_table)

    # Save and download ANOVA results
    anova_results_file = "ANOVAResults.csv"
    anova_table.to_csv(anova_results_file)
    st.download_button(
        label="📥 Download ANOVA Results",
        data=open(anova_results_file, "rb").read(),
        file_name="ANOVAResults.csv",
        mime="text/csv"
    )

    # Shapiro-Wilk Test
    stat, p_shapiro = shapiro(data["Gene"])
    st.write(f"📊 **Shapiro-Wilk Test**: p-value = `{p_shapiro:.5f}`")

    # Kruskal-Wallis Test
    groups = [data[data["structure_color"] == region]["Gene"] for region in data["structure_color"].unique()]
    stat, p_kruskal = kruskal(*groups)
    st.write(f"📊 **Kruskal-Wallis Test**: p-value = `{p_kruskal:.5f}`")

    # Dunn's Post Hoc Test (if Kruskal-Wallis is significant)
    if p_kruskal < 0.05:
        st.success("✅ The Kruskal-Wallis test is significant! Performing Dunn’s post hoc test...")
        dunn_result = sp.posthoc_dunn(data, val_col="Gene", group_col="structure_color", p_adjust="bonferroni")

        st.write("### 📊 Dunn’s Post Hoc Test Results")
        st.dataframe(dunn_result)

        # Save and download Dunn's test results
        dunn_results_file = "DunnResults.csv"
        dunn_result.to_csv(dunn_results_file)
        st.download_button(
            label="📥 Download Dunn's Test Results",
            data=open(dunn_results_file, "rb").read(),
            file_name="DunnResults.csv",
            mime="text/csv"
        )
    else:
        st.info("ℹ️ Kruskal-Wallis test is not significant. No post hoc test needed.")

