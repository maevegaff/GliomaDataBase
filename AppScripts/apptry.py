import streamlit as st
import pandas as pd
from statsmodels.formula.api import ols
import statsmodels.api as sm
from scipy.stats import shapiro, kruskal
import scikit_posthocs as sp
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy.stats as stats
from scipy.stats import shapiro
from scipy.stats import kruskal
import scikit_posthocs as sp

def process_file(uploaded_file):
    # Load the file uploaded by the user
    data = pd.read_csv(uploaded_file)

    # Display the data
    st.write(data.head())

    # Summary statistics
    summary_stats = data.describe()
    st.write(summary_stats)

    # Save the summary statistics to a CSV file
    summary_stats_file = 'SumStat.csv'
    summary_stats.to_csv(summary_stats_file)
    create_download_button(summary_stats_file, "Download Summary Statistics")

    # Define the model
    model = ols('Gene ~ C(structure_color)', data=data).fit()

    # Perform ANOVA
    anova_table = sm.stats.anova_lm(model, typ=2)
    st.write(anova_table)

    # Save the ANOVA results to a CSV file
    anova_results_file = 'ANOVAResults.csv'
    anova_table.to_csv(anova_results_file)
    create_download_button(anova_results_file, "Download ANOVA Results")

    # Shapiro-Wilk Test
    gene_column = 'Gene'  # Change this to the actual column name if different
    stat, p = shapiro(data[gene_column])
    st.write(f"Shapiro-Wilk Test: p-value = {p}")

    # Run Kruskal-Wallis Test
    groups = [data[data['structure_color'] == region][gene_column] for region in data['structure_color'].unique()]
    stat, p = kruskal(*groups)
    st.write(f"Kruskal-Wallis Test: p-value = {p}")

    # If significant, run Dunn’s test
    if p < 0.05:
        dunn = sp.posthoc_dunn(data, val_col=gene_column, group_col='structure_color', p_adjust='bonferroni')
        st.write(dunn)

        # Save the results to a CSV file
        dunn_test_results_file = dunn
        dunn.to_csv(dunn_test_results_file)
        create_download_button(dunn_test_results_file, "Download Dunn's Test Results")
    else:
        st.write("Dunn's test was not significant, so no results file was created.")

def create_download_button(file_path, label):
    with open(file_path, "rb") as file:
        st.download_button(
            label=label,
            data=file.read(),
            file_name=file_path,
            mime="text/csv"
        )

st.title('Tumor Region Expression Analysis: Statistical Tests and Results')

uploaded_file = st.file_uploader("Choose a CSV file", type="csv")

if uploaded_file is not None:
    # Process the file
    process_file(uploaded_file)

