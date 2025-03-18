import streamlit as st
import pandas as pd
import os
import glob
from statsmodels.formula.api import ols
import statsmodels.api as sm
from scipy.stats import shapiro, kruskal
import scikit_posthocs as sp
import seaborn as sns
import matplotlib.pyplot as plt

def process_files(file1_path):
    # Load the file
    file1 = pd.read_csv(file1_path)

    # Save the summary statistics to a CSV file
    summary_stats = file1.describe()
    summary_stats.to_csv('SumStat.csv')

    # Define the model
    model = ols('Gene ~ C(structure_color)', data=file1).fit()

    # Perform ANOVA
    anova_table = sm.stats.anova_lm(model, typ=2)

    # Print the results
    print(anova_table)

    # Save the ANOVA results to a CSV file
    anova_table.to_csv('ANOVAResults.csv')

    # Shapiro-Wilk Test
    stat, p = shapiro(file1['Gene'])  # Replace 'expression' with actual column
    print(f"Shapiro-Wilk Test: p-value = {p}")


    # Run Kruskal-Wallis Test
    stat, p = kruskal(*[file1[file1['structure_color'] == region]['Gene'] for region in file1['structure_color'].unique()])
    print(f"Kruskal-Wallis Test: p = {p}")

    # If significant, run Dunn’s test
    dunn = sp.posthoc_dunn(file1, val_col='Gene', group_col='structure_color', p_adjust='bonferroni')
    print(dunn)

    # Save the results to a CSV file
    dunn.to_csv('DunnResults.csv')

  

    # Perform Dunn's post hoc test if Kruskal-Wallis test is significant
    if p < 0.05:
        st.write('Dunn\'s post hoc test result:')
        st.write(dunn)

        # Create a box plot
        plt.figure(figsize=(10, 6))
        sns.boxplot(x='Gene', y='structure_color', data=file1)

        # Precompute max and min values for each gene
        gene_max = file1.groupby('Gene')['structure_color'].max()
        gene_min = file1.groupby('Gene')['structure_color'].min()

        # Mark significance on the plot
        significant_pairs = dunn[dunn < 0.05].stack().index.tolist()
        for (gene1, gene2) in significant_pairs:
            y_max = max(gene_max[gene1], gene_max[gene2])
            y_min = min(gene_min[gene1], gene_min[gene2])
            plt.plot([gene1, gene2], [y_max + 0.1, y_max + 0.1], 'k-', lw=1.5)
            plt.text((gene1 + gene2) / 2, y_max + 0.15, '*', ha='center', va='bottom', color='k')

        plt.title('Box plot with significant differences marked')
        plt.xlabel('Gene')
        plt.ylabel('Structure Color')
        plt.xticks(rotation=45)
        plt.tight_layout()

        # Save the plot
        plot_path = 'box_plot.png'
        plt.savefig(plot_path)
        st.write('Box plot with significant differences:')
        st.image(plot_path)
    else:
        st.write('Kruskal-Wallis test is not significant, skipping Dunn\'s post hoc test.')

    return plot_path

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

if uploaded_file1 is not None:
    # Save the uploaded file to disk
    file1_path = os.path.join("uploaded_file1.csv")
    
    with open(file1_path, "wb") as f:
        f.write(uploaded_file1.getbuffer())
    
    # Process the file
    new_file_path = process_files(file1_path)
    
    # Display the link to download the new file
    st.success(f"File processed successfully! [Download box_plot.png](./{new_file_path})")
    st.download_button(
        label="Download box_plot.png",
        data=open(new_file_path, "rb").read(),
        file_name="box_plot.png",
        mime="image/png"
    )
