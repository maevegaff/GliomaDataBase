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
    file2 = pd.read_csv(uploaded_file2, header=None, index_col=0, encoding='utf-8')  # Read first row as data

    # Debug: Check the contents of the loaded files
    st.write("### File 1 Preview")
    st.dataframe(file1.head())
    st.write("### File 2 Preview")
    st.dataframe(file2.head())

    # Check if file2 is empty
    if file2.empty:
        st.error("❌ The Gene Expression CSV file is empty.")
        st.stop()

    # Transpose and reshape file2 for correct merging
    file2_transposed = file2.T.reset_index()
    file2_transposed.rename(columns={'index': 'donor_id'}, inplace=True)
    
    # Debug: Check the contents of the transposed file2
    st.write("### Transposed File 2 Preview")
    st.dataframe(file2_transposed.head())

    # Check if file2_transposed is empty
    if file2_transposed.empty:
        st.error("❌ The transposed Gene Expression DataFrame is empty.")
        st.stop()

    file2_long = file2_transposed.melt(id_vars=['donor_id'], var_name='Gene', value_name='Expression')
    
    # Debug: Check the contents of the transposed and reshaped file2
    st.write("### Transposed and Reshaped File 2 Preview")
    st.dataframe(file2_long.head())
    
    # Convert donor_id to string for merging
    file2_long['donor_id'] = file2_long['donor_id'].astype(str)
    file1['donor_id'] = file1['donor_id'].astype(str)
    
    # Select required columns from file1
    required_columns = ['donor_id', 'structure_color', 'survival_days']
    if not all(col in file1.columns for col in required_columns):
        st.error(f"❌ The IvyTumorInformation CSV must contain these columns: {required_columns}")
        st.stop()
    
    # Merge datasets on donor_id
    merged_df = pd.merge(file1[required_columns], file2_long, on='donor_id', how='inner')
    
    # Debug: Check the contents of the merged DataFrame
    st.write("### Merged DataFrame Preview")
    st.dataframe(merged_df.head())
    
    # Check if merged_df is empty
    if merged_df.empty:
        st.error("❌ The merged DataFrame is empty. Please check the input files for matching donor IDs.")
        st.stop()
    
    # Filter only rows where `structure_color == '05D004'`
    merged_df = merged_df[merged_df['structure_color'] == '05D004']
    merged_df.drop(columns=['structure_color'], inplace=True)
    
    # Debug: Check the contents of the filtered DataFrame
    st.write("### Filtered DataFrame Preview")
    st.dataframe(merged_df.head())
    
    # Convert survival_days to numeric and remove missing values
    merged_df['survival_days'] = pd.to_numeric(merged_df['survival_days'], errors='coerce')
    merged_df.dropna(subset=['survival_days'], inplace=True)
    
    # Debug: Check the contents of the DataFrame after removing missing values
    st.write("### DataFrame After Removing Missing Values")
    st.dataframe(merged_df.head())
    
    # Set event column to 1 (assuming death occurs at survival_days)
    merged_df['event'] = np.where(merged_df['survival_days'] > 0, 1, 0)
    
    # Compute mean gene expression per donor
    mean_expression_per_donor = merged_df.groupby('donor_id', as_index=False)['Expression'].mean()
    merged_df = pd.merge(merged_df, mean_expression_per_donor, on='donor_id', suffixes=('', '_Mean'))
    
    # Compute median expression for High/Low classification
    median_expression = merged_df['Expression_Mean'].median()
    merged_df['Gene_High_Low'] = np.where(merged_df['Expression_Mean'] >= median_expression, 'High Expression', 'Low Expression')
    
    # Display merged dataset preview
    st.write("### 🔍 Merged Data Preview")
    st.dataframe(merged_df.head())
    
    # Save the merged file
    merged_file_path = 'merged_file.csv'
    merged_df.to_csv(merged_file_path, index=False, float_format='%.4f')
    st.success(f"Merged file saved as {merged_file_path}")
    
    # Download processed data
    csv_data = merged_df.to_csv(index=False).encode('utf-8')
    st.download_button("📥 Download Processed Data", data=csv_data, file_name="Processed_Survival_Data.csv", mime="text/csv")
    
    # Fit Cox Model for survival analysis
    df = merged_df[['survival_days', 'event', 'Gene_High_Low']]
    df['Gene_High_Low'] = df['Gene_High_Low'].astype('category').map({'Low Expression': 0, 'High Expression': 1})
    
    if df['survival_days'].nunique() > 1 and df['Gene_High_Low'].nunique() > 1:
        cph = CoxPHFitter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                cph.fit(df, duration_col='survival_days', event_col='event', formula="Gene_High_Low")
                st.write("The Cox Proportional Hazards Model summary provides statistical details about the relationship between gene expression levels and survival time.")
            except Exception as e:
                st.error(f"❌ An error occurred while fitting the Cox model: {e}")
                st.stop()
        
        st.write("### 📊 Cox Proportional Hazards Model Summary")
        st.text(cph.print_summary())
    else:
        st.warning("⚠️ Not enough variability in the data to fit the Cox model.")
    
    # Kaplan-Meier Survival Curve
    st.write("### 📈 Kaplan-Meier Survival Curve for High vs. Low Gene Expression")
    kmf_low = KaplanMeierFitter()
    kmf_high = KaplanMeierFitter()
    
    low_expression = merged_df[merged_df['Gene_High_Low'] == 'Low Expression']
    high_expression = merged_df[merged_df['Gene_High_Low'] == 'High Expression']
    
    if len(low_expression) > 0 and len(high_expression) > 0:
        try:
            kmf_low.fit(low_expression['survival_days'], event_observed=low_expression['event'], label="Low Expression")
        except Exception as e:
            st.error(f"❌ An error occurred while fitting the Kaplan-Meier model for Low Expression: {e}")
            st.stop()
        
        try:
            kmf_high.fit(high_expression['survival_days'], event_observed=high_expression['event'], label="High Expression")
        except Exception as e:
            st.error(f"❌ An error occurred while fitting the Kaplan-Meier model for High Expression: {e}")
            st.stop()
    
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
        st.download_button("📥 Download Kaplan-Meier Plot", data=buf, file_name="Kaplan_Meier_Survival_Plot.png", mime="image/png")
    else:
        st.warning("⚠️ One of the groups has no survival data, so the Kaplan-Meier plot may not display properly.")
