import streamlit as st
import pandas as pd
import numpy as np
import io
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
import warnings

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
    file2 = pd.read_csv(uploaded_file2, index_col=0)  # Assume first column is an index

    # Transpose file2 (Genes as columns, donor IDs as rows)
    file2_transposed = file2.T.reset_index()
    file2_transposed.rename(columns={'index': 'Gene'}, inplace=True)
    
    # Convert donor_id to string to ensure proper merging
    file2_transposed['Gene'] = file2_transposed['Gene'].astype(str)
    file1['donor_id'] = file1['donor_id'].astype(str)
    
    # Select required columns from file1
    required_columns = ['donor_id', 'structure_color', 'survival_days']
    if not all(col in file1.columns for col in required_columns):
        st.error(f"❌ The IvyTumorInformation CSV must contain these columns: {required_columns}")
        st.stop()

    # Merge datasets
    merged_df = pd.merge(file1[required_columns], file2_transposed, on='Gene', how='inner')
    
    # Filter only rows where `structure_color == '05D004'`
    merged_df = merged_df[merged_df['structure_color'] == '05D004']
    
    # Drop structure_color column as it's no longer needed
    merged_df.drop(columns=['structure_color'], inplace=True)
    
    # Convert survival_days to numeric
    merged_df['survival_days'] = pd.to_numeric(merged_df['survival_days'], errors='coerce')
    
    # Remove rows with missing survival days
    merged_df.dropna(subset=['survival_days'], inplace=True)
    
    # Set event column to 1 (since death always occurs after survival days)
    merged_df['event'] = 1
    
    # Compute mean gene expression per donor
    gene_columns = [col for col in merged_df.columns if col not in ['donor_id', 'survival_days', 'event']]
    mean_expression_per_donor = merged_df.groupby('donor_id', as_index=False)[gene_columns].mean()
    
    # Merge mean expression back into main dataframe
    merged_df = pd.merge(merged_df[['donor_id', 'survival_days', 'event']], mean_expression_per_donor, on='donor_id', how='left')
    
    # Compute median gene expression for one representative gene (modify as needed)
    representative_gene = gene_columns[0]  # Use the first gene as an example
    median_expression = merged_df[representative_gene].median()
    
    # Create High/Low Expression groups
    merged_df['Gene_High_Low'] = np.where(merged_df[representative_gene] >= median_expression, 'High Expression', 'Low Expression')
    
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
    df['Gene_High_Low'] = df['Gene_High_Low'].astype('category').map({'Low Expression': 0, 'High Expression': 1})
    
    # Ensure there is variability in survival times
    if df['survival_days'].nunique() > 1 and df['Gene_High_Low'].nunique() > 1:
        cph = CoxPHFitter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cph.fit(df, duration_col='survival_days', event_col='event', formula="Gene_High_Low")
        
        # Display Cox Model Summary
        st.write("### 📊 Cox Proportional Hazards Model Summary")
        st.text(cph.print_summary())
    else:
        st.warning("⚠️ Not enough variability in the data to fit the Cox model.")
    
    # Kaplan-Meier Survival Curve
    st.write("### 📈 Kaplan-Meier Survival Curve for High vs. Low Gene Expression")
    
    kmf_low = KaplanMeierFitter()
    kmf_high = KaplanMeierFitter()
    
    # Split data into Low and High Expression groups
    low_expression = merged_df[merged_df['Gene_High_Low'] == 'Low Expression']
    high_expression = merged_df[merged_df['Gene_High_Low'] == 'High Expression']
    
    if len(low_expression) > 0 and len(high_expression) > 0:
        kmf_low.fit(low_expression['survival_days'], event_observed=low_expression['event'], label="Low Expression")
        kmf_high.fit(high_expression['survival_days'], event_observed=high_expression['event'], label="High Expression")

        fig, ax = plt.subplots(figsize=(8, 6))
        kmf_low.plot(ax=ax, ci_show=True)
        kmf_high.plot(ax=ax, ci_show=True)

        ax.set_title('Kaplan-Meier Survival Curve for High vs. Low Gene Expression')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Survival Probability')
        ax.legend()

        st.pyplot(fig)
        
        buf = io.BytesIO()
        fig.savefig(buf, format='png')
        buf.seek(0)
        st.download_button(
            label="📥 Download Kaplan-Meier Plot",
            data=buf,
            file_name="Kaplan_Meier_Survival_Plot.png",
            mime="image/png"
        )
    else:
        st.warning("⚠️ One of the groups has no survival data, so the Kaplan-Meier plot may not display properly.")
