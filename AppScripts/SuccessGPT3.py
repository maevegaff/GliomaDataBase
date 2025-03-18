
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
    columns_to_select = ['donor_id', 'donor_name', 'structure_color', 'tumor_name', 'molecular_subtype']
    file1_selected = file1[columns_to_select]

    # Merge files
    merged_df = pd.concat([file1_selected, file2_transposed], axis=1)

    # Save processed data
    processed_file = "ProcessedData.csv"
    merged_df.to_csv(processed_file, index=False, float_format="%.4f")
    
    return processed_file, merged_df

# Function to generate and save boxplot with significance markers
def generate_boxplot(data):
    plt.figure(figsize=(10, 6))
    ax = sns.boxplot(x="structure_color", y="Gene", data=data, palette="Set2")
    
    # Perform Dunn's post hoc test
    dunn_result = sp.posthoc_dunn(data, val_col="Gene", group_col="structure_color", p_adjust="bonferroni")

    # Get unique tumor regions
    regions = data["structure_color"].unique()
    
    # Define mapping from color codes to readable labels
    color_to_label = {
        '218FA5': 'Leading Edge',
        'D104D0': 'Infiltrating Tumour',
        '05D004': 'Cellular Tumor',
        '43D1F8': 'Perinecrotic Zone',
        '05D0AA': 'Pseudopalisading',
        'FF6600': 'HP Blood Vessel',
        'FF3300': 'MV Proliferation'
    }
    
    # Change x-axis labels
    ax.set_xticklabels([color_to_label.get(label, label) for label in regions])


    # Replace structure_color codes with readable names
    data['Region'] = data['structure_color'].map(color_to_label).fillna(data['structure_color'])
    
    # Define y-position for significance markers
    y_max = data["Gene"].max()  
    y_offset = (y_max - data["Gene"].min()) * 0.05  
    y_pos = y_max + y_offset  

    # Iterate over all region pairs and annotate significance
    for (i, j) in itertools.combinations(range(len(regions)), 2):
        p_value = dunn_result.iloc[i, j]  
        if p_value < 0.05:  # Significant
            x1, x2 = i, j  
            ax.plot([x1, x1, x2, x2], [y_pos, y_pos + y_offset, y_pos + y_offset, y_pos], color="black", linewidth=1)
            ax.text((x1 + x2) / 2, y_pos + y_offset * 1.2, "*", ha="center", va="bottom", fontsize=14, color="red")
            y_pos += y_offset * 1.5  

    plt.xticks(rotation=45)
    plt.xlabel("Tumor Region")
    plt.ylabel("Gene Expression")
    plt.title("Gene Expression Across Tumor Regions")

    # Save and display the plot
    plot_filename = "boxplot.png"
    plt.savefig(plot_filename, bbox_inches="tight")
    st.pyplot(plt)  # Display in Streamlit
    
    return plot_filename

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

    # Ensure required columns exist
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

    # Generate and display boxplot
    st.write("### 📊 Gene Expression Boxplot")
    boxplot_file = generate_boxplot(data)

    # Download button for boxplot
    with open(boxplot_file, "rb") as f:
        st.download_button(
            label="📥 Download Boxplot",
            data=f,
            file_name="GeneExpressionBoxplot.png",
            mime="image/png"
        )


