import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

def load_and_process_files(file1_path, file2_path):
    file1 = pd.read_csv(file1_path)
    file2 = pd.read_csv(file2_path)
    file2 = file2.iloc[:, 1:].transpose().reset_index()
    file2.columns = ['Gene'] + [f'Sample_{i}' for i in range(1, len(file2.columns))]
    file2['Gene'] = file2['Gene'].astype(str)
    file1['donor_id'] = file1['donor_id'].astype(str)
    merged_df = pd.merge(file1, file2, left_on='donor_id', right_on='Gene', how='inner')
    return merged_df

def run_kaplan_meier_analysis(data):
    data['survival_days'] = pd.to_numeric(data['survival_days'], errors='coerce')
    data = data.dropna(subset=['survival_days', 'mgmt_methylation'])
    data['mgmt_methylation'] = data['mgmt_methylation'].astype(str)
    kmf = KaplanMeierFitter()
    plt.figure(figsize=(10, 6))
    for group in data['mgmt_methylation'].unique():
        mask = data['mgmt_methylation'] == group
        kmf.fit(data['survival_days'][mask], event_observed=np.ones(sum(mask)))
        kmf.plot(label=f'MGMT Methylation: {group}')
    plt.title('Kaplan-Meier Survival Curve')
    plt.xlabel('Survival Days')
    plt.ylabel('Survival Probability')
    st.pyplot(plt)

st.title('Kaplan-Meier Survival Analysis')
file1 = st.file_uploader("Upload Patient Data (CSV)", type="csv")
file2 = st.file_uploader("Upload Gene Expression Data (CSV)", type="csv")
if file1 and file2:
    with open("file1.csv", "wb") as f:
        f.write(file1.getbuffer())
    with open("file2.csv", "wb") as f:
        f.write(file2.getbuffer())
    data = load_and_process_files("file1.csv", "file2.csv")
    st.write("### Merged Data Preview")
    st.dataframe(data.head())
    run_kaplan_meier_analysis(data)
