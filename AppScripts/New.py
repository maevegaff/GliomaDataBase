import streamlit as st
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import io
import matplotlib.pyplot as plt

# Streamlit UI
st.title('🧬 Gene Expression & Survival Analysis')

st.markdown("""
The app will:
- Merge the datasets.
- Classify samples into **High** and **Low Expression** groups.
- Fit a **Cox Proportional Hazards Model**.
- Plot **survival curves** for both groups.
""")

# Upload files
uploaded_file1 = st.file_uploader("📂 Upload IvyTumorInformation CSV", type="csv")
uploaded_file2 = st.file_uploader("📂 Upload Gene Expression CSV", type="csv")

if uploaded_file1 and uploaded_file2:
    # Load files
    file1 = pd.read_csv(uploaded_file1)
    file2 = pd.read_csv(uploaded_file2)

    # Remove first column from file2 (assumed to be index-like)
    file2 = file2.iloc[:, 1:]

    # Transpose file2
    file2_transposed = file2.transpose().reset_index()
    file2_transposed.columns = ['Gene'] + [f'Sample_{i}' for i in range(1, len(file2_transposed.columns))]

    # Ensure numeric values
    for col in file2_transposed.columns[1:]:
        file2_transposed[col] = pd.to_numeric(file2_transposed[col], errors='coerce')

    # Select relevant columns from file1
    required_columns = ['donor_id', 'donor_name', 'structure_color', 'tumor_name', 'molecular_subtype', 'survival_days']
    if not all(col in file1.columns for col in required_columns):
        st.error(f"❌ The IvyTumorInformation CSV must contain these columns: {required_columns}")
        st.stop()

    # Merge datasets
    merged_df = pd.concat([file1[required_columns], file2_transposed], axis=1)

    # Remove rows without numeric survival days
    merged_df = merged_df[pd.to_numeric(merged_df['survival_days'], errors='coerce').notnull()]

    # Filter rows where structure_color is not equal to '05D004'
    merged_df = merged_df[merged_df['structure_color'] == '05D004']

    # Ensure Gene column is numeric
    merged_df['Gene'] = pd.to_numeric(merged_df['Gene'], errors='coerce')

    # Calculate median gene expression for each donor_id using vectorized operations
    median_expression_per_donor = merged_df.groupby('donor_id')['Gene'].median().reset_index()
    median_expression_per_donor.columns = ['donor_id', 'Median_Gene_Expression']

    # Merge median expression back into the main dataframe
    merged_df = pd.merge(merged_df, median_expression_per_donor, on='donor_id', how='left')

    # Update Gene column with median expression for each donor_id
    merged_df['Gene'] = merged_df['Median_Gene_Expression']

    # Add event column
    # Assumption: All patients had the event (e.g., death) for the purpose of this analysis.
    merged_df['event'] = 1  

    # Ensure Gene column is numeric
    merged_df['Gene'] = pd.to_numeric(merged_df['Gene'], errors='coerce')

    # Compute median gene expression
    median_expression = merged_df['Gene'].median()

    # Categorize into High and Low Expression groups
    merged_df['Gene_Group'] = np.where(merged_df['Gene'] >= median_expression, 'High Expression', 'Low Expression')

    # Display merged dataset preview
    st.write("### 🔍 Merged Data Preview")
    st.dataframe(merged_df.head())

    # Allow downloading processed data
    csv_data = merged_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Download Processed Data",
        data=csv_data,
        file_name="Processed_Survival_Data.csv",
        mime="text/csv"
    )

    # Prepare data for Cox Proportional Hazards model
    df = merged_df[['survival_days', 'Gene_Group', 'event']]
    df['Gene_Group'] = df['Gene_Group'].map({'Low Expression': 0, 'High Expression': 1})

    # Fit Cox Model
    cph = CoxPHFitter()
    cph.fit(df, duration_col='survival_days', event_col='event', formula="Gene_Group")

    # Display Cox Model Summary
    st.write("### 📊 Cox Proportional Hazards Model Summary")
    st.text(cph.print_summary())

    # Plot survival function
    st.write("### 📈 Survival Function for High vs. Low Gene Expression")

    # Create the survival plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Calculate and plot survival function for high and low expression groups
    for group in ['Low Expression', 'High Expression']:
        temp_df = df[df['Gene_Group'] == (1 if group == 'High Expression' else 0)]
        cph.plot_partial_effects_on_outcome(
            covariates='Gene_Group', 
            values=[temp_df['Gene_Group'].mean()], 
            ax=ax, 
            label=group
        )

    # Customize plot
    ax.set_title('Survival Function for High vs. Low Gene Expression')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Survival Probability')
    ax.legend()
    
    # Show plot in Streamlit
    st.pyplot(fig)

    # Allow downloading the plot
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    st.download_button(
        label="📥 Download Survival Plot",
        data=buf,
        file_name="Survival_Plot.png",
        mime="image/png"
    )

    # Display statistics for each gene
    st.write("### 📊 Gene Expression Statistics")

    # Calculate statistics
    stats = file2.describe().transpose()
    stats['IQR'] = stats['75%'] - stats['25%']
    stats = stats[['min', '25%', '50%', '75%', 'max', 'IQR']]
    stats.columns = ['Min', '25th Percentile', 'Median', '75th Percentile', 'Max', 'IQR']

    # Display statistics
    st.dataframe(stats)

    # Allow downloading statistics
    stats_csv = stats.to_csv(index=True).encode('utf-8')
    st.download_button(
        label="📥 Download Gene Expression Statistics",
        data=stats_csv,
        file_name="Gene_Expression_Statistics.csv",
        mime="text/csv"
    )