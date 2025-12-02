import os
import pandas as pd

# File path for the samples list
samples_file = './data/samples.txt'

# Base directory for the experiment files
BASE_DIR = '/data/metaEdit/experiments/eggnog'
base_dir = f"{BASE_DIR}/RNA_edit_result/RNA_filtered/"

eggnog_dir = f"{BASE_DIR}/eggnog_out/"

# Read SAMPLE_IDs from the samples.txt file
samples = [line.strip() for line in open(samples_file, 'r') if line.strip()]

# Iterate through each SAMPLE_ID
for sample_id in samples:
    rna_file_path = os.path.join(base_dir, sample_id, f"{sample_id}_filtered_annotated.tsv")
    eggnog_file_path = os.path.join(eggnog_dir, sample_id, ".emapper.annotations")

    # Check if the files exist
    if not os.path.isfile(rna_file_path):
        print(f"RNA file not found for SAMPLE_ID: {sample_id} at {rna_file_path}")
        continue

    if not os.path.isfile(eggnog_file_path):
        print(f"EggNOG file not found for SAMPLE_ID: {sample_id} at {eggnog_file_path}")
        continue

    try:
        # Read the RNA file
        rna_df = pd.read_csv(rna_file_path, sep='\t')

        # Extract the first gene ID from the 'Gene_id' column
        rna_df['First_Gene_id'] = rna_df['Gene_id'].str.split(',').str[0]

        # Read the EggNOG file
        eggnog_df = pd.read_csv(eggnog_file_path, sep='\t', low_memory=False, skiprows=4)

        # Merge RNA and EggNOG dataframes on the first gene ID and #query
        merged_df = rna_df.merge(eggnog_df, left_on='First_Gene_id', right_on='#query', how='left')

        # Output the merged dataframe or save it
        output_path = os.path.join(base_dir, sample_id, f"{sample_id}_merged_output.tsv")
        merged_df.to_csv(output_path, sep='\t', index=False)

        print(f"Merged file saved for SAMPLE_ID: {sample_id} at {output_path}")

    except Exception as e:
        print(f"Error processing files for SAMPLE_ID: {sample_id}. Error: {e}")
