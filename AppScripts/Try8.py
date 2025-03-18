import streamlit as st
import subprocess
import os

def run_analysis(file_path):
    try:
    st.text_area("Standard Output", result.stdout, height=200)
    st.text_area("Standard Error", result.stderr, height=200)
        st.text(result.stderr)
    except subprocess.CalledProcessError as e:
        st.error(f"An error occurred: {e}")

st.title('Tumor Region Stats Analysis')
st.write('Upload a CSV file and click the button below to run the analysis.')

uploaded_file = st.file_uploader("Choose a CSV file", type="csv")

if uploaded_file is not None:
    file_path = os.path.join("/tmp", uploaded_file.name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    st.write(f"File uploaded: {uploaded_file.name}")

    if st.button('Run Analysis'):
        run_analysis(file_path)