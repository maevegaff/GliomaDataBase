import streamlit as st
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt

# Streamlit UI
st.title('🧬 Gene Expression & Survival Analysis')

st.markdown("""
Upload a CSV file containing **gene expression** and **survival days** data.  
The app will:
- Classify samples into **High** and **Low Expression** groups based on the median.
- Fit a **Cox Proportional Hazards Model**.
- Plot **survival curves** for both groups.
""")

# Upload file
uploaded_file = st.file_uploader("📂 Upload your CSV file", type="csv")

if uploaded_file:
    # Load Data
    data = pd.read_csv(uploaded_file)

    # Ensure necessary columns exist
    if 'survival_days' not in data.columns or 'Gene' not in data.columns:
        st.error("❌ The CSV file must contain 'survival_days' and 'Gene' columns.")
        st.stop()

    # Remove rows without valid survival days
    data = data[pd.to_numeric(data['survival_days'], errors='coerce').notnull()]

    # Add 'event' column assuming all patients experienced the event
    data['event'] = 1  

    # Compute the median gene expression
    median_expression = data['Gene'].median()

    # Categorize into High and Low Expression groups
    data['Gene_Group'] = np.where(data['Gene'] >= median_expression, 'High Expression', 'Low Expression')

    # Display dataset preview
    st.write("### 🔍 Data Preview")
    st.dataframe(data.head())

    # Allow downloading processed data
    csv_data = data.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Download Processed Data",
        data=csv_data,
        file_name="Processed_Survival_Data.csv",
        mime="text/csv"
    )

    # Prepare data for Cox Proportional Hazards model
    df = data[['survival_days', 'Gene_Group', 'event']]
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
