import streamlit as st
import pandas as pd
from scipy.stats import kruskal
import scikit_posthocs as sp

st.title("Test Streamlit App")
st.write("If you see this message, Streamlit is working!")

uploaded_file = st.file_uploader("Upload a CSV file", type="csv")

if uploaded_file is not None:
    st.success("File uploaded successfully!")
    st.write("Filename:", uploaded_file.name)

    # Read the CSV file into a DataFrame
    df = pd.read_csv(uploaded_file)
    
    # Display the DataFrame
    st.write(df)
    
    # Perform Kruskal-Wallis H-test on 'structure_color' and 'Gene' columns
    if 'structure_color' in df.columns and 'Gene' in df.columns:
        groups = [df[df['structure_color'] == g]['Gene'].dropna().values for g in df['structure_color'].unique()]
        stat, p = kruskal(*groups)
        
        st.write(f"Kruskal-Wallis H-test statistic: {stat}")
        st.write(f"P-value: {p}")
        
        # Perform Dunn's test if Kruskal-Wallis test is significant
        if p < 0.05:
            st.write("Performing Dunn's test...")
            dunn_results = sp.posthoc_dunn(df, val_col='Gene', group_col='structure_color', p_adjust='bonferroni')
            st.write(dunn_results)
            
            # Save Dunn's test results to a CSV file
            dunn_results_file = "dunn_results.csv"
            dunn_results.to_csv(dunn_results_file)
            st.success(f"Dunn's test results saved to {dunn_results_file}")
            
            # Provide a download button for the results file
            with open(dunn_results_file, "rb") as file:
                st.download_button(
                    label="Download Dunn's test results",
                    data=file,
                    file_name=dunn_results_file,
                    mime="text/csv"
                )
        else:
            st.write("Kruskal-Wallis test is not significant, skipping Dunn's test.")
    else:
        st.write("The required columns 'structure_color' and 'Gene' are not present in the uploaded file.")