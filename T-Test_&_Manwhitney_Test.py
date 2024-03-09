#Script will check  test for normality (shapiro-wilk test) and and based on the distribution chhose the statistical test
#(t-test or Mann Whitney U test)
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind, shapiro
import openpyxl
# Read data from CSV file
df = pd.read_csv("top10_gene_count.txt", sep='\t')

# List of genes of Interest
genes = ["MUC3A", "MUC16", "TTN", "OBSCN", "FSIP2", "CCDC168", "LOC400499", "CCDC187", "DNAH14", "MK67", "CDH1"]

# Perform statistical test for each gene
results = {}
for gene in genes:
    # Check normality
    shapiro_gc = shapiro(df[df["Condition"] == "Gastric Cancer"][gene])
    shapiro_hc = shapiro(df[df["Condition"] == "Healthy Control"][gene])

    gc_data = df[df["Condition"] == "Gastric Cancer"][gene]
    hc_data = df[df["Condition"] == "Healthy Control"][gene]

    # Calculate means for each group
    gc_mean = gc_data.mean()
    hc_mean = hc_data.mean()

    # If data is normally distributed, perform t-test
    if shapiro_gc.pvalue > 0.05 and shapiro_hc.pvalue > 0.05:
        t_stat, t_p_value = ttest_ind(df[df["Condition"] == "Gastric Cancer"][gene],
                                      df[df["Condition"] == "Healthy Control"][gene])
        if t_p_value < 0.05:
            p_value_significance = str(t_p_value) + "*"
        else:
            p_value_significance = str(t_p_value)
        results[gene] = {"Statistic": t_stat, "p-value": p_value_significance, "Test": "t", "gc_mean_(n=13)": gc_mean, "hc_mean_(n=26)": hc_mean}
    # If data is not normally distributed, perform Mann-Whitney U test
    else:
        mwu_stat, mwu_p_value = mannwhitneyu(df[df["Condition"] == "Gastric Cancer"][gene],
                                             df[df["Condition"] == "Healthy Control"][gene])
        if mwu_p_value < 0.05:
            p_value_significance = str(mwu_p_value) + "*"
        else:
            p_value_significance = str(mwu_p_value)
        results[gene] = {"Statistic": mwu_stat, "p-value": p_value_significance, "Test": "MWU", "gc_mean_(n=13)": gc_mean, "hc_mean_(n=26)": hc_mean}

# Display results
results_df = pd.DataFrame.from_dict(results, orient='index')
print(results_df)
print("\n\n\n\n\n")
results_df.to_excel("gene_analysis_results.xlsx")
