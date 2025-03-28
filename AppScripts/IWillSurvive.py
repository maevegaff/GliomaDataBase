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
    file2_transposed = file2.T.reset_index(drop=True)
    
    # Debug: Check the contents of the transposed file2
    st.write("### Transposed File 2 Preview")
    st.dataframe(file2_transposed.head())

    # Check if file2_transposed is empty
    if file2_transposed.empty:
        st.error("❌ The transposed Gene Expression DataFrame is empty.")
        st.stop()

    # Flatten the transposed DataFrame to a single column
    file2_long = file2_transposed.melt(var_name='Gene', value_name='Expression')
    
    # Debug: Check the contents of the transposed and reshaped file2
    st.write("### Transposed and Reshaped File 2 Preview")
    st.dataframe(file2_long.head())
    
    # Add the Expression column to file1
    file1['Expression'] = file2_long['Expression'].values
    
    # Debug: Check the contents of the updated file1
    st.write("### Updated File 1 with Expression Column")
    st.dataframe(file1.head())
    
    # Select required columns from file1
    required_columns = ['donor_name', 'structure_color', 'survival_days', 'Expression']
    if not all(col in file1.columns for col in required_columns):
        st.error(f"❌ The IvyTumorInformation CSV must contain these columns: {required_columns}")
        st.stop()
    
    # Filter only rows where `structure_color == '05D004'`
    filtered_df = file1[file1['structure_color'] == '05D004']
    filtered_df.drop(columns=['structure_color'], inplace=True)
    
    # Debug: Check the contents of the filtered DataFrame
    st.write("### Filtered DataFrame Preview")
    st.dataframe(filtered_df.head())
    
    # Convert survival_days to numeric and remove missing values
    filtered_df['survival_days'] = pd.to_numeric(filtered_df['survival_days'], errors='coerce')
    filtered_df.dropna(subset=['survival_days'], inplace=True)
    
    # Debug: Check the contents of the DataFrame after removing missing values
    st.write("### DataFrame After Removing Missing Values")
    st.dataframe(filtered_df.head())
    
    # Set event column to 1 (assuming death occurs at survival_days)
    filtered_df['event'] = np.where(filtered_df['survival_days'] > 0, 1, 0)
    
    # Compute mean gene expression per donor
    mean_expression_per_donor = filtered_df.groupby('donor_name', as_index=False)['Expression'].mean()
    filtered_df = pd.merge(filtered_df, mean_expression_per_donor, on='donor_name', suffixes=('', '_Mean'))
    
    # Keep the first row in each group and filter out the rest
    filtered_df = filtered_df.groupby('donor_name').first().reset_index()
    
    # Debug: Check the contents of the DataFrame after keeping the first row in each group
    st.write("### DataFrame After Keeping First Row in Each Group")
    st.dataframe(filtered_df.head())

    # Compute median expression for High/Low classification
    median_expression = filtered_df['Expression_Mean'].median()
    filtered_df['Gene_High_Low'] = np.where(filtered_df['Expression_Mean'] >= median_expression, 'High Expression', 'Low Expression')
    
    # Display merged dataset preview
    st.write("### 🔍 Merged Data Preview")
    st.dataframe(filtered_df.head())
    
    # Save the merged file
    merged_file_path = 'merged_file.csv'
    filtered_df.to_csv(merged_file_path, index=False, float_format='%.4f')
    st.success(f"Merged file saved as {merged_file_path}")
    
    # Download processed data
    csv_data = filtered_df.to_csv(index=False).encode('utf-8')
    st.download_button("📥 Download Processed Data", data=csv_data, file_name="Processed_Survival_Data.csv", mime="text/csv")
    
    # Fit Cox Model for survival analysis
    df = filtered_df[['survival_days', 'event', 'Gene_High_Low']]
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
    
    low_expression = filtered_df[filtered_df['Gene_High_Low'] == 'Low Expression']
    high_expression = filtered_df[filtered_df['Gene_High_Low'] == 'High Expression']
    
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
        ax.set_title('Kaplan-Meier Survival Curve for High vs. Low 5HT3A Expression')
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
