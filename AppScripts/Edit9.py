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
    file1 = pd.read_csv(uploaded_file1, encoding='utf-8')
    file2 = pd.read_csv(uploaded_file2, index_col=0, encoding='utf-8')  # Read first column as index

    # Debug: Check file contents
    st.write("### File 1 Preview")
    st.dataframe(file1.head())
    st.write("### File 2 Preview")
    st.dataframe(file2.head())
    
    # Check if file2 is empty
    if file2.empty:
        st.error("❌ The Gene Expression CSV file is empty.")
        st.stop()
    
    # Transpose file2 to align donors properly
    file2_transposed = file2.T.reset_index()
    file2_transposed.rename(columns={'index': 'donor_id'}, inplace=True)
    
    # Debug: Check transposed file
    st.write("### Transposed File 2 Preview")
    st.dataframe(file2_transposed.head())
    
    # Ensure donor_id is a string for merging
    file1['donor_id'] = file1['donor_name'].astype(str)
    file2_transposed['donor_id'] = file2_transposed['donor_id'].astype(str)
    
    # Reshape file2 for merging
    file2_long = file2_transposed.melt(id_vars=['donor_id'], var_name='Gene', value_name='Expression')
    
    # Debug: Check reshaped file2
    st.write("### Reshaped Gene Expression Data Preview")
    st.dataframe(file2_long.head())
    
    # Validate required columns in file1
    required_columns = ['donor_id', 'survival_days']
    if not all(col in file1.columns for col in required_columns):
        st.error(f"❌ The IvyTumorInformation CSV must contain these columns: {required_columns}")
        st.stop()
    
    # Merge datasets
    merged_df = pd.merge(file1[['donor_id', 'survival_days']], file2_long, on='donor_id', how='inner')
    
    # Debug: Check merged DataFrame
    st.write("### Merged DataFrame Preview")
    st.dataframe(merged_df.head())
    
    if merged_df.empty:
        st.error("❌ No matching donor IDs found in both datasets.")
        st.stop()
    
    # Convert survival_days to numeric and clean data
    merged_df['survival_days'] = pd.to_numeric(merged_df['survival_days'], errors='coerce')
    merged_df.dropna(subset=['survival_days'], inplace=True)
    
    # Set event column to 1 (assuming death occurs at survival_days)
    merged_df['event'] = np.where(merged_df['survival_days'] > 0, 1, 0)
    
    # Compute mean gene expression per donor
    mean_expression = merged_df.groupby('donor_id', as_index=False)['Expression'].mean()
    merged_df = pd.merge(merged_df, mean_expression, on='donor_id', suffixes=('', '_Mean'))
    
    # Classify gene expression into High/Low groups
    median_expression = merged_df['Expression_Mean'].median()
    merged_df['Gene_High_Low'] = np.where(merged_df['Expression_Mean'] >= median_expression, 'High Expression', 'Low Expression')
    
    # Save merged dataset
    csv_data = merged_df.to_csv(index=False).encode('utf-8')
    st.download_button("📥 Download Processed Data", data=csv_data, file_name="Processed_Survival_Data.csv", mime="text/csv")
    
    # Cox Proportional Hazards Model
    df = merged_df[['survival_days', 'event', 'Gene_High_Low']]
    df['Gene_High_Low'] = df['Gene_High_Low'].astype('category').cat.codes
    
    if df['survival_days'].nunique() > 1 and df['Gene_High_Low'].nunique() > 1:
        cph = CoxPHFitter()
        try:
            cph.fit(df, duration_col='survival_days', event_col='event', formula="Gene_High_Low")
            st.write("### 📊 Cox Proportional Hazards Model Summary")
            st.text(cph.print_summary())
        except Exception as e:
            st.error(f"❌ Error fitting Cox model: {e}")
    else:
        st.warning("⚠️ Not enough variability to fit the Cox model.")
    
    # Kaplan-Meier Survival Curve
    st.write("### 📈 Kaplan-Meier Survival Curve")
    kmf_low, kmf_high = KaplanMeierFitter(), KaplanMeierFitter()
    
    low_expr, high_expr = merged_df[merged_df['Gene_High_Low'] == 'Low Expression'], merged_df[merged_df['Gene_High_Low'] == 'High Expression']
    
    fig, ax = plt.subplots(figsize=(8, 6))
    if not low_expr.empty:
        kmf_low.fit(low_expr['survival_days'], event_observed=low_expr['event'], label="Low Expression")
        kmf_low.plot(ax=ax, ci_show=True)
    if not high_expr.empty:
        kmf_high.fit(high_expr['survival_days'], event_observed=high_expr['event'], label="High Expression")
        kmf_high.plot(ax=ax, ci_show=True)
    
    ax.set_title('Kaplan-Meier Survival Curve')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Survival Probability')
    ax.legend()
    
    st.pyplot(fig)
    
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    st.download_button("📥 Download Kaplan-Meier Plot", data=buf, file_name="Kaplan_Meier_Survival_Plot.png", mime="image/png")
