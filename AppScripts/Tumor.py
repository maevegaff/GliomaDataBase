# Streamlit interface
st.title('Glioma Data Analysis')

# File uploader
uploaded_file = st.file_uploader("Choose a CSV file", type="csv")

if uploaded_file is not None:
    # Load the data
    data = pd.read_csv(uploaded_file)
    
    # Display the data
    st.write("Data Preview:")
    st.write(data.head())
    
    # Summary statistics
    st.write("Summary Statistics:")
    summary_stats = data.describe()
    st.write(summary_stats)
    
    # Save the summary statistics to a CSV file
    summary_stats.to_csv('SumStat.csv')
    
    # Define the model
    model = ols('Gene ~ C(structure_color)', data=data).fit()
    
    # Perform ANOVA
    anova_table = sm.stats.anova_lm(model, typ=2)
    
    # Display the results
    st.write("ANOVA Results:")
    st.write(anova_table)
    
    # Save the ANOVA results to a CSV file
    anova_table.to_csv('ANOVAResults.csv')
    
    # Shapiro-Wilk Test
    stat, p = shapiro(data['Gene'])
    st.write(f"Shapiro-Wilk Test: p-value = {p}")
    
    # Run Kruskal-Wallis Test
    stat, p = kruskal(*[data[data['structure_color'] == region]['Gene'] for region in data['structure_color'].unique()])
    st.write(f"Kruskal-Wallis Test: p = {p}")
    
    # If significant, run Dunn’s test
    if p < 0.05:
        dunn = sp.posthoc_dunn(data, val_col='Gene', group_col='structure_color', p_adjust='bonferroni')
        st.write("Dunn's Test Results:")
        st.write(dunn)
        
        # Save the results to a CSV file
        dunn.to_csv('DunnResults.csv')
