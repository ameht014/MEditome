# Extract orthologs of interest
import pandas as pd

BASE_DIR = '/data/metaEdit'
orth_file = f"{BASE_DIR}/experiments/MSA_merge/MSA_merged/MSA_summary_aggregated_filter5.tsv"

df = pd.read_csv(orth_file, sep='\t')
df = df[['seed_ortholog','alignment_position']].drop_duplicates()
df.to_csv('data/extracted_orthologs_filter5.csv', index=False)

