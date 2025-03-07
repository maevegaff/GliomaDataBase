import pandas as pd
import seaborn as sns
from scipy.stats import kruskal
from scikit_posthocs import posthoc_dunn
import matplotlib.pyplot as plt

# Load the data file
data_file = 'new_file.csv'
data = pd.read_csv(data_file)

# Load the Dunn test results
dunn_results = pd.read_csv('DunnResults.csv')

# Perform Kruskal-Wallis test
kruskal_results = kruskal(*[data[data['structure_color'] == region]['Gene'] for region in data['structure_color'].unique()])

# Perform Dunn's post hoc test
dunn_results = posthoc_dunn(data, val_col='Gene', group_col='structure_color', p_adjust='bonferroni')



# Create a box plot comparing gene expression across 7 different tumour regions
plt.figure(figsize=(12, 6))
sns.boxplot(x='structure_color', y='Gene', data=data)
plt.title('Gene Expression Across Different Tumour Regions')
plt.xlabel('Tumour Region')
plt.ylabel('HTR1A Gene Expression')

# Change the labels on the x-axis
ax = plt.gca()
labels = [item.get_text() for item in ax.get_xticklabels()]
labels = ['Leading Edge' if label == '218FA5' else
          'Infiltrating Tumour' if label == 'D104D0' else
          'Cellular Tumor' if label == '05D004' else
          'Perinecrotic Zone' if label == '43D1F8' else
          'Pseudopalisading' if label == '05D0AA' else
          'HP Blood Vessel' if label == 'FF6600' else
          'MV Proliferation' if label == 'FF3300' else label
          for label in labels]
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=6)


# Add significance markers to the plot
def add_significance(ax, dunn_results, data):
    box_pairs = [(region1, region2) for region1 in data['structure_color'].unique() for region2 in data['structure_color'].unique() if region1 < region2]
    y_max = data['Gene'].max()
    y_offset = y_max * 0.15
    stagger_offset = y_offset * 0.75  # Increased stagger offset for better visibility
    for i, (region1, region2) in enumerate(box_pairs):
        p_value = dunn_results.loc[region1, region2]
        if p_value < 0.05:
            x1, x2 = data['structure_color'].unique().tolist().index(region1), data['structure_color'].unique().tolist().index(region2)
            y = y_max + y_offset + (i % 2) * stagger_offset
            h = y_offset
            ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, color='k')
            ax.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color='k')

ax = sns.boxplot(x='structure_color', y='Gene', data=data)
add_significance(ax, dunn_results, data)
plt.show()