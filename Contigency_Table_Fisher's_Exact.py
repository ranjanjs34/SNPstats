import os
import csv
from tqdm import tqdm
from scipy.stats import fisher_exact

def create_contingency_table(variant_file, healthy_directory, cancer_directory):
    # Read the variant list from the file
    variants = []
    with open(variant_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')  # Assuming tab-separated variant file
            variant = '_'.join(columns[:5])  # Combine the first three columns with underscores as the variant name
            variants.append(variant)

    # Initialize the contingency table with headers
    contingency_table = [[''] + variants, ['Healthy'] + [0] * len(variants), ['Gastric Cancer'] + [0] * len(variants)]

    # Match and count occurrences in each healthy data file
    for filename in tqdm(os.listdir(healthy_directory), desc='Processing Healthy Data'):
        filepath = os.path.join(healthy_directory, filename)
        with open(filepath, 'r') as file:
            for line in file:
                data = line.strip().split('\t')  # Assuming tab-separated data file
                variant = '_'.join(data[:5])  # Combine the first three columns with underscores as the variant name
                if variant in variants:
                    idx = variants.index(variant) + 1  # Get the index of the variant in the contingency table
                    contingency_table[1][idx] += 1  # Increment count for healthy samples with the variant

    # Match and count occurrences in each gastric cancer data file
    for filename in tqdm(os.listdir(cancer_directory), desc='Processing Gastric Cancer Data'):
        filepath = os.path.join(cancer_directory, filename)
        with open(filepath, 'r') as file:
            for line in file:
                data = line.strip().split('\t')  # Assuming tab-separated data file
                variant = '_'.join(data[:5])  # Combine the first three columns with underscores as the variant name
                if variant in variants:
                    idx = variants.index(variant) + 1  # Get the index of the variant in the contingency table
                    contingency_table[2][idx] += 1  # Increment count for cancer samples with the variant

    return contingency_table

def save_contingency_table(contingency_table, output_file, delimiter=','):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=delimiter)
        writer.writerows(contingency_table)

def perform_fisher_exact_test(contingency_table):
    p_values = []
    for i in range(1, len(contingency_table[0])):
        oddsratio, p_value = fisher_exact([[contingency_table[1][i], contingency_table[2][i]],
                                           [sum(contingency_table[1][1:]), sum(contingency_table[2][1:])]])
        p_values.append(p_value)
    return p_values
# Example usage
variant_file = 'variant_list.txt'  # Path to the file containing the variant list
healthy_directory = '/mnt/c/Users/Ranjan/Desktop/All_GC/filtered/GC'  # Path to the directory containing the healthy data files
cancer_directory = '/mnt/c/Users/Ranjan/Desktop/All_GC/filtered/Healthy'  # Path to the directory containing the gastric cancer data files
output_file = 'contingency_table.csv'  # Path to the output file for the contingency table
p_value_file = 'p_values.csv'
# Create the contingency table
contingency_table = create_contingency_table(variant_file, healthy_directory, cancer_directory)

# Save the contingency table as CSV file
save_contingency_table(contingency_table, output_file, delimiter=',')

print("Contingency table saved as", output_file)
p_values = perform_fisher_exact_test(contingency_table)

# Save p-values as tab-delimited text file

with open(p_value_file, 'w') as file:
    header = '\t'.join(['Variant', "Fisher's Exact p-value"])
    file.write(f"{header}\n")
    for variant, p_value in zip(contingency_table[0][1:], p_values):
        file.write(f"{variant}\t{p_value}\n")

print("P-values saved as", p_value_file)
