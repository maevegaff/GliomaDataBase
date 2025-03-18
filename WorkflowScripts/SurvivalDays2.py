import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv(r'C:\Users\maeve\.vscode\GliomaDataBase\new_file.csv')

# Ensure the data has the necessary columns
if 'survival_days' not in data.columns or 'Gene' not in data.columns:
    raise ValueError("The CSV file must contain 'survival_days' and 'Gene' columns.")
    
# Remove rows without a numeric value in 'survival_days' column
data = data[pd.to_numeric(data['survival_days'], errors='coerce').notnull()]

# Add the event column (assuming all patients experienced the event)
data['event'] = 1  

# Compute the median gene expression
median_expression = data['Gene'].median()

# Categorize gene expression into "High" and "Low"
data['Gene_Group'] = np.where(data['Gene'] >= median_expression, 'High Expression', 'Low Expression')

# Prepare the data for Cox Proportional Hazards model
df = data[['survival_days', 'Gene_Group', 'event']]

# Convert categorical variable to numerical encoding
df['Gene_Group'] = df['Gene_Group'].map({'Low Expression': 0, 'High Expression': 1})

# Fit the Cox model
cph = CoxPHFitter()
cph.fit(df, duration_col='survival_days', event_col='event', formula="Gene_Group")

# Print model summary
cph.print_summary()

# Plot survival function for high and low expression groups
plt.figure(figsize=(8, 6))

# Calculate survival function for high and low expression
for group in ['Low Expression', 'High Expression']:
    temp_df = df[df['Gene_Group'] == (1 if group == 'High Expression' else 0)]
    cph.plot_partial_effects_on_outcome(covariates='Gene_Group', values=[temp_df['Gene_Group'].mean()], label=group)

# Customize plot
plt.title('Survival Function for High and Low Gene Expression')
plt.xlabel('Time (days)')
plt.ylabel('Survival Probability')
plt.legend()
plt.show()
