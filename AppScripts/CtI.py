import streamlit as st
import pandas as pd
from scipy.stats import kruskal
import scikit_posthocs as sp

def main():
    st.title("Tumor Region Expression Analysis: Statistical Tests and Results")

    uploaded_file = st.file_uploader("Choose a CSV file", type="csv")

    if uploaded_file is not None:
        data = pd.read_csv(uploaded_file)

        if data.shape[1] != 2:
            st.error("CSV file must have exactly two columns")
            return

        groups = data.iloc[:, 0]
        values = data.iloc[:, 1]

        kruskal_result = kruskal(*[values[groups == g] for g in groups.unique()])
        st.write(f"Kruskal-Wallis Test: p-value = {kruskal_result.pvalue}")

        if kruskal_result.pvalue < 0.05:
            dunn_result = sp.posthoc_dunn(data, val_col=data.columns[1], group_col=data.columns[0], p_adjust='bonferroni')
            st.write("Kruskal-Wallis test is significant. Dunn Post Hoc Test results:")
            st.dataframe(dunn_result)

            # Save the Dunn's test results to a CSV file
            dunn_results_file = 'DunnResults.csv'
            dunn_result.to_csv(dunn_results_file)
            with open(dunn_results_file, "rb") as file:
                st.download_button(
                    label="Download Dunn's Test Results",
                    data=file.read(),
                    file_name=dunn_results_file,
                    mime="text/csv"
                )
        else:
            st.info("Kruskal-Wallis test is not significant, no need for post hoc test")

if __name__ == "__main__":
    main()
