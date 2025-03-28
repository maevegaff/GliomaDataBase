import streamlit as st
import pandas as pd
import numpy as np
import io
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter

# Streamlit UI
st.title('🧬 Gene Expression & Survival Analysis')

st.markdown("""
The app will:
- Merge the datasets.
- Filter for structure color **05D004**.
- Average gene expression values for duplicate donor IDs.
- Fit a **Cox Proportional Hazards Model**.
- Plot a **Kaplan-Meier survival curve** for high vs. low gene expression groups.
""")

# Upload files
uploaded_file1 = st.file_uploader("📂 Upload IvyTumorInformation CSV", type="csv")
uploaded_file2 = st.file_uploader("📂 Upload Gene Expression CSV", type="csv")

if uploaded_file1 and uploaded_file2:
    # Load files
    file1 = pd.read_csv(uploaded_file1)
    file2 = pd.read_csv(uploaded_file2)

    # Remove first column from file2 (assuming it's an index column)
    file2 = file2.iloc[:, 1:]

    # Transpose file2
    file2_transposed = file2.transpose().reset_index()
    file2_transposed.columns = ['Gene'] + [f'Sample_{i}' for i in range(1, len(file2_transposed.columns))]

    # Convert columns to numeric
    for col in file2_transposed.columns[1:]:
        file2_transposed[col] = pd.to_numeric(file2_transposed[col], errors='coerce')

    # Select relevant columns from file1
    required_columns = ['donor_id', 'donor_name', 'structure_color', 'tumor_name', 'molecular_subtype', 'survival_days']
    if not all(col in file1.columns for col in required_columns):
        st.error(f"❌ The IvyTumorInformation CSV must contain these columns: {required_columns}")
        st.stop()

    # Merge datasets
    merged_df = pd.concat([file1[required_columns], file2_transposed], axis=1)

    # Remove rows with non-numeric survival days
    merged_df = merged_df[pd.to_numeric(merged_df['survival_days'], errors='coerce').notnull()]

    # Filter only rows where `structure_color == '05D004'`
    merged_df = merged_df[merged_df['structure_color'] == '05D004']

    # Convert 'Gene' column to numeric
    merged_df['Gene'] = pd.to_numeric(merged_df['Gene'], errors='coerce')

    # Compute mean gene expression per donor_id
    mean_expression_per_donor = merged_df.groupby('donor_id', as_index=False)['Gene'].mean()

    # Merge mean expression back into the main dataframe
    merged_df = merged_df.drop(columns=['Gene'])  # Drop the original column
    merged_df = pd.merge(merged_df, mean_expression_per_donor, on='donor_id', how='left')

    # Drop duplicate donor_id rows (keep only first occurrence)
    merged_df = merged_df.drop_duplicates(subset=['donor_id'], keep='first')

    # Compute median gene expression
    median_expression = merged_df['Gene'].median()

    # Create a new column for High/Low Expression groups
    merged_df['Gene_High_Low'] = np.where(merged_df['Gene'] >= median_expression, 'High Expression', 'Low Expression')

    # Create event column: 1 if above median, 0 if below
    merged_df['event'] = np.where(merged_df['Gene'] >= median_expression, 1, 0)

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

    # Fit Cox Model for statistical analysis
    df = merged_df[['survival_days', 'event', 'Gene_High_Low']]
    df['Gene_High_Low'] = df['Gene_High_Low'].map({'Low Expression': 0, 'High Expression': 1})

    cph = CoxPHFitter()
    cph.fit(df, duration_col='survival_days', event_col='event', formula="Gene_High_Low")

    # Display Cox Model Summary
    st.write("### 📊 Cox Proportional Hazards Model Summary")
    st.text(cph.print_summary())

    # Create Kaplan-Meier survival curves
    st.write("### 📈 Kaplan-Meier Survival Curve for High vs. Low Gene Expression")

    kmf_low = KaplanMeierFitter()
    kmf_high = KaplanMeierFitter()

    # Split data into Low and High Expression groups
    low_expression = merged_df[merged_df['Gene_High_Low'] == 'Low Expression']
    high_expression = merged_df[merged_df['Gene_High_Low'] == 'High Expression']

    # Fit Kaplan-Meier estimators
    kmf_low.fit(low_expression['survival_days'], event_observed=low_expression['event'], label="Low Expression")
    kmf_high.fit(high_expression['survival_days'], event_observed=high_expression['event'], label="High Expression")

    # Plot Kaplan-Meier survival curves
    fig, ax = plt.subplots(figsize=(8, 6))
    kmf_low.plot(ax=ax, ci_show=False)
    kmf_high.plot(ax=ax, ci_show=False)

    # Customize plot
    ax.set_title('Kaplan-Meier Survival Curve for High vs. Low Gene Expression')
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
        label="📥 Download Kaplan-Meier Plot",
        data=buf,
        file_name="Kaplan_Meier_Survival_Plot.png",
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
