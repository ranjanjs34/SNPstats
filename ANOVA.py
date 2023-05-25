import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

# Read the data from the CSV file
data = pd.read_csv('MUC3A_NonSyno_Count.csv')

# Separate the data into cancer and healthy groups
cancer_samples = data[data['Condition'] == 'Gastric Cancer']['Number_of_MUC3A_Missence_Variants']
healthy_samples = data[data['Condition'] == 'Healthy Control']['Number_of_MUC3A_Missence_Variants']

# Calculate mean, standard deviation, and other statistics
cancer_mean = cancer_samples.mean()
healthy_mean = healthy_samples.mean()

cancer_std = cancer_samples.std()
healthy_std = healthy_samples.std()

# Perform One-Way ANOVA
statistic, p_value = stats.f_oneway(cancer_samples, healthy_samples)
corrected_p_value = multipletests(p_value, method='fdr_bh')[1]

# Create a DataFrame with the results
results_df = pd.DataFrame({
    'Group': ['Gastric Cancer', 'Healthy Control'],
    'Mean': [cancer_mean, healthy_mean],
    'Standard Deviation': [cancer_std, healthy_std],
    'F-statistic': [statistic, ''],
    'p-value': [p_value, ''],
    'FDR corrected p-value': [corrected_p_value, '']
})

# Save the DataFrame as a tab-delimited CSV file
results_df.to_csv('output_file.csv', sep='\t', index=False)
