import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv(r'C:\Users\maeve\.vscode\GliomaDataBase\new_file.csv')

# Print the CSV file
print(data)

# Ensure the data has the necessary columns
if 'survival_days' not in data.columns or 'Gene' not in data.columns:
    raise ValueError("The CSV file must contain 'survival_days' and 'Gene' columns.")
    
# Remove rows without a numeric value in 'survival_days' column
data = data[pd.to_numeric(data['survival_days'], errors='coerce').notnull()]

 # Print the cleaned dataframe
print(data)

# Add the event column with all values set to 1
data['event'] = 1

# Prepare the data for the Cox Proportional Hazards model
df = data[['survival_days', 'Gene', 'event']]

# Fit the Cox Proportional Hazards model
cph = CoxPHFitter()
cph.fit(df, duration_col='survival_days', event_col='event', formula="Gene")

# Print the summary of the model
cph.print_summary()



# Plot survival function for different gene expression levels
cph.plot_partial_effects_on_outcome(covariates='Gene', values=np.linspace(data['Gene'].min(), data['Gene'].max(), 100))



# Show the plot
plt.title('Survival Function for Different Gene Expression Levels')
plt.xlabel('Time (days)')
plt.ylabel('Survival Probability')
plt.show()
