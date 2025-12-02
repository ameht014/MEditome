#!/usr/bin/env python3

import os
from tqdm import tqdm
import pandas as pd
from multiprocessing import Pool
import argparse
import time

## We are looking for A-I RNA editing.
## Therefore, the change should be either A to G (because Inosine is interpreted as guanosine by the molecular machineries)
## Or


def get_edit_level_annotated(all_position_info, result_file, pos_list, cov_threshold, edit_threshold):
    out_edits_dict = dict()
    out_cov_dict = dict()
    out_prot_dict = dict()

    if os.path.exists(result_file) and os.path.exists(all_position_info)\
            and os.path.getsize(result_file)>0 and os.path.getsize(all_position_info)>0:

        ## Step 1 - merge all_position_info and _result.txt
        position_info_df = pd.read_csv(all_position_info, sep='\t', names=["Accession", "Position", "Base",
                                                                             "cov", "covA", "covC", "covG", "covT"])

        result_df = pd.read_csv(result_file, sep='\t')

        merged_df = result_df.merge(position_info_df, on=["Position","Accession"], how='left')
        merged_df['Position'] = merged_df['Position'].apply(lambda x: str(x))
        merged_df['Position'] = merged_df['Accession'] + '_' + merged_df['Position']
        merged_df = merged_df[merged_df['Position'].isin(pos_list)]

        # ## Iterate through each position and report edit level, coverage, and annotation

        for index, row in merged_df.iterrows():
            pos = row['Position']
            covA, covC, covG, covT = row['covA'], row['covC'], row['covG'], row['covT']
            refBase = row['Base']
            edit_level = get_edit_level_one(refBase, covA, covC, covG, covT)
            if int(row['cov']) > cov_threshold and edit_level>edit_threshold:
                out_edits_dict[pos] = edit_level
            out_cov_dict[pos] = (covA, covC, covG, covT)
            out_prot_dict[pos] = (row['Old_base'], row['New_base'], row['Gene_biotype'], row['Gene_id'])
    return out_edits_dict, out_cov_dict, out_prot_dict



def get_edit_level_one(refBase, covA, covC, covG, covT):
    if (refBase == "A" or refBase == "a") and float(covG + covA)>0:
        # Disregard edits A to T and A to C
        edit_level = int(covG) / float(covG + covA)

    elif (refBase == "T" or refBase == "t") and float(covC + covT)>0:
        # Disregard edits T to G and T to A
        edit_level = int(covC) / float(covC + covT)

    else: # If reference not A or T
        edit_level = 0

    return  edit_level

# def get_edit_level(all_position_info, cov_threshold=30, edit_threshold=0, pos_list=None, get_cov=True):
#     out_edits_dict = dict()
#     out_cov_dict = dict()
#     if os.path.exists(all_position_info):
#         with open(all_position_info, 'r') as f:
#             for line in f.readlines():
#                 pass
#                 id, pos, refBase, cov, covA, covC, covG, covT = line.strip('\n').split('\t')
#                 if pos_list is not None and pos not in pos_list:
#                     pass
#                 #pdb.set_trace()
#                 if pos_list is not None and pos not in pos_list:
#                     continue
#                 if get_cov:
#                     out_cov_dict[pos] = (covA, covC, covG, covT)
#                     pass
#
#                 if int(cov)>cov_threshold:
#                     edit_level = get_edit_level_one(refBase, covA, covC, covG, covT)
#                     # if pos == '3607424' and int(cov) > cov_threshold:
#                     #     print(edit_level)
#                     if edit_level > float(edit_threshold):
#                         out_edits_dict[pos] = edit_level
#                         pass
#
#     return out_edits_dict, out_cov_dict


def get_edit_level(all_position_info, cov_threshold=30, edit_threshold=0, pos_list=None, get_cov=True):
    # Ensure the file exists
    if not os.path.exists(all_position_info):
        print(f"{all_position_info} DOESN't")
        return {}, {}  # Return empty dicts if file does not exist

    # Read the file into a DataFrame
    df = pd.read_csv(all_position_info, sep='\t', header=None,
                     names=['id', 'pos', 'refBase', 'cov', 'covA', 'covC', 'covG', 'covT'])
    df['pos'] =  df['id'] + '_' + df['pos'].astype(str)
    # Filter by pos_list if it is not None
    if pos_list is not None:
        df = df[df['pos'].isin([str(pos) for pos in pos_list])]

    # Initialize dictionaries
    out_edits_dict = {}
    out_cov_dict = {}

    # Calculate coverage if get_cov is True and convert pos to string explicitly if not already
    if get_cov:
        out_cov_dict = df.set_index('pos')[['covA', 'covC', 'covG', 'covT']].T.to_dict('list')

    # Assuming get_edit_level_one is defined and can work with Series or individual values
    if 'get_edit_level_one' in globals():
        # Filter rows where coverage is above threshold and make a copy to avoid SettingWithCopyWarning
        high_cov_df = df[df['cov'] > cov_threshold].copy()

        if len(high_cov_df)>0:
            # Calculate edit levels for rows with high coverage
            high_cov_df['edit_level'] = high_cov_df.apply(
                lambda row: get_edit_level_one(row['refBase'], row['covA'], row['covC'], row['covG'], row['covT']),
                axis=1
            )

            # Filter by edit threshold
            filtered_df = high_cov_df[high_cov_df['edit_level'] > edit_threshold]
            out_edits_dict = filtered_df.set_index('pos')['edit_level'].to_dict()

    return out_edits_dict, out_cov_dict

# def write_edit_to_file(unique_positions,all_dict, out_file, all_samples):
#     #unique_positions = list(unique_positions)
#     print(f"{len(unique_positions)} unique positions.")
#     print(f"Writing to a file {out_file}")
#     with open(out_file, 'w') as out:
#         out.write('sample' + ',' + ','.join(unique_positions)+'\n')
#         for i, sample in enumerate(all_samples):
#             next_line = []
#             next_line.append(sample)
#             for pos in unique_positions:
#                 if pos in all_dict[i].keys():
#                     next_line.append(str(all_dict[i][pos]))
#                 else:
#                     next_line.append(str(0))
#             out.write(','.join(next_line)+'\n')

def write_edit_to_file(unique_positions, all_dict, out_file, all_samples):
    print(f"{len(unique_positions)} unique positions.")
    print(f"Writing to a file {out_file}")

    # Prepare data for transposed output
    # Initialize a dictionary to hold values for each position across all samples
    transposed_data = {str(pos): [] for pos in unique_positions}
    for pos in unique_positions:
        pos = str(pos)
        for i in range(len(all_samples)):
            if pos in all_dict[i].keys():
                transposed_data[pos].append(str(all_dict[i][pos]))
            else:
                transposed_data[pos].append('0')

    # Write to the file
    with open(out_file, 'w') as out:
        # Write the header row with sample names
        out.write('position,' + ','.join(all_samples) + '\n')
        # Write each row for the positions
        for pos, values in transposed_data.items():
            out.write(pos + ',' + ','.join(values) + '\n')


# def write_cov_to_file(unique_positions, cov_dict, out_file, all_samples):
#     #unique_positions = list(unique_positions)
#     print(f"{len(unique_positions)} unique positions.")
#     print(f"Writing to a file {out_file}")
#     all_cols = []
#     # covA, covC, covG, covT
#     for pos in unique_positions:
#         all_cols.append(pos+"_A")
#         all_cols.append(pos + "_C")
#         all_cols.append(pos + "_G")
#         all_cols.append(pos + "_T")
#
#     with open(out_file, 'w') as out:
#         out.write('sample' + ',' + ','.join(all_cols) + '\n')
#         for i, sample in enumerate(all_samples):
#             next_line = []
#             next_line.append(sample)
#             for pos in unique_positions:
#                 if pos in cov_dict[i].keys():
#                     for cov in cov_dict[i][pos]:
#                         next_line.append(str(cov))
#                 else:
#                     # Write 0 to A,T,C, and G
#                     next_line.append(str(0))
#                     next_line.append(str(0))
#                     next_line.append(str(0))
#                     next_line.append(str(0))
#             out.write(','.join(next_line) + '\n')

def write_cov_to_file(unique_positions, cov_dict, out_file, all_samples):
    print(f"{len(unique_positions)} unique positions.")
    print(f"Writing to a file {out_file}")

    # Prepare header row with sample names, starting with a placeholder for the position column
    header_row = ['position'] + all_samples

    # Prepare data for transposed output
    transposed_rows = []
    for pos in unique_positions:
        for nt in ['A', 'C', 'G', 'T']:  # Nucleotides
            row_key = f"{pos}_{nt}"
            row = [row_key]  # Start each row with the position_nucleotide key
            for i in range(len(all_samples)):
                if pos in cov_dict[i].keys():
                    # Assuming cov_dict[i][pos] is a list with 4 elements: [covA, covC, covG, covT]
                    nt_index = ['A', 'C', 'G', 'T'].index(nt)
                    row.append(str(cov_dict[i][pos][nt_index]))
                else:
                    row.append('0')  # If position not in sample, append '0'
            transposed_rows.append(row)

    # Write to the file
    with open(out_file, 'w') as out:
        out.write(','.join(header_row) + '\n')  # Write the header
        for row in transposed_rows:
            out.write(','.join(row) + '\n')  # Write each transposed row

    print(f"Transposed coverage data written to {out_file}.")


def write_annot_to_file(all_dict_prot_annot, out_rna_annot):
    rows = []
    for pos in all_dict_prot_annot.keys():
        rows.append([pos] + list(all_dict_prot_annot[pos]))
    df = pd.DataFrame(rows, columns=['Position', 'Old_base', 'New_base', 'Gene_biotype', 'Gene_id'])
    df.to_csv(out_rna_annot, sep='\t', index=False)
    return 1

# def compute_coverage_sample(sample, rna_dir, ref_name, rna_cov_threshold, edit_threshold, unique_positions):
#     out_rna_dir = rna_dir + '/' + sample + '.' + ref_name + '/'
#     all_position_info_rna = out_rna_dir + sample + '.' + ref_name + '_all_position_info.txt'
#     result_file = out_rna_dir + sample + '.' + ref_name + '_result.txt'
#     if not os.path.exists(result_file):
#         out_rna_dir = rna_dir + '/'
#         all_position_info_rna = out_rna_dir + sample + '_all_position_info.txt'
#         result_file = out_rna_dir + sample + '_result.txt'
#

def process_sample_dna(args):
    sample, dna_dir, ref_name, cov_threshold, edit_threshold, all_rna_sites = args
    out_dna_dir = dna_dir + '/' + sample + '.' + ref_name + '/'
    all_position_info = out_dna_dir + sample + '.' + ref_name + '_all_position_info.txt'
    if not os.path.exists(all_position_info):
        out_dna_dir = dna_dir + '/'
        all_position_info = out_dna_dir + sample + '_all_position_info.txt'
    edits_dict, cov_dict = get_edit_level(all_position_info, cov_threshold, edit_threshold, pos_list=all_rna_sites)
    return edits_dict, cov_dict

def get_dna_snps_whole_study(dna_dir, samples_list, ref_name, rna_edit_study, out_snps, out_cov, cov_threshold=30, edit_threshold=0.03, n_cores=32, out_dir='./'):
    """
    Get dna_edit dictionary across the whole study.
        Keep the largest edit score.
        This will be used for finale RNA editing step
    The script will create a file for a given ref_name where columns are positions, rows are samples, and values
                                                                    are percentages of variant bases (edit_score)

    """
    # all_dna_folders = os.listdir(dna_folder)
    # for dna_folder_i in all_dna_folders:
    # result_file = os.path.join(dna_folder, dna_folder, dna_folder + "_result.txt")
    #all_samples = [x.strip('\n') for x in open(samples_list).readlines()]

    all_samples = samples_list
    print("Here")
    print(rna_edit_study)
    all_rna_df = pd.read_csv(rna_edit_study)
    all_rna_sites = list(all_rna_df['position'].apply(lambda x: str(x)).unique())
    #all_rna_df=None
    # with open(rna_edit_study, 'r') as f:
    #     all_rna_sites = f.readline().strip('\n').split(',')[1:]

    print(f"Reading DNA edit from {len(all_samples)} samples", flush=True)
    print(f"Checking {len(all_rna_sites)} potential sites", flush=True)
    # args = [(sample, dna_dir, ref_name, cov_threshold, edit_threshold, all_rna_sites) for sample in all_samples]
    # with Pool(processes=n_cores) as pool:
    #     results = list(tqdm(pool.imap(process_sample_dna, args), total=len(all_samples)))

    # ## Sequential:
    # out_dir = out_dir + '/sites/'
    # if os.path.exists(out_dir):
    #     shutil.rmtree(out_dir)
    # os.mkdir(out_dir)

    results = []
    for sample in tqdm(all_samples):
        result = process_sample_dna(args=(sample, dna_dir, ref_name, cov_threshold, edit_threshold, all_rna_sites))
        #pdb.set_trace()
        results.append(result)
    print(f"Done. Writing to files...", flush=True)

    all_dict_snps = [result[0] for result in results]
    all_dict_cov = [result[1] for result in results]

    # Write to file
    write_edit_to_file(all_rna_sites, all_dict_snps, out_snps, all_samples)
    write_cov_to_file(all_rna_sites, all_dict_cov, out_cov, all_samples)
    return 1


def compute_editing_sample(sample, rna_dir, ref_name, rna_cov_threshold, edit_threshold):
    out_rna_dir = rna_dir + '/' + sample + '.' + ref_name + '/'
    all_position_info_rna = out_rna_dir + sample + '.' + ref_name + '_all_position_info.txt'
    if not os.path.exists(all_position_info_rna):
        out_rna_dir = rna_dir + '/'
        all_position_info_rna = out_rna_dir + sample + '_all_position_info.txt'
    #print(all_position_info_rna)
    rna_edit_dict, _ = get_edit_level(all_position_info_rna, rna_cov_threshold, edit_threshold, get_cov=False)
    return rna_edit_dict

def compute_coverage_sample(sample, rna_dir, ref_name, rna_cov_threshold, edit_threshold, unique_positions):
    out_rna_dir = rna_dir + '/' + sample + '.' + ref_name + '/'
    all_position_info_rna = out_rna_dir + sample + '.' + ref_name + '_all_position_info.txt'
    result_file = out_rna_dir + sample + '.' + ref_name + '_result.txt'
    if not os.path.exists(result_file):
        out_rna_dir = rna_dir + '/'
        all_position_info_rna = out_rna_dir + sample + '_all_position_info.txt'
        result_file = out_rna_dir + sample + '_result.txt'

    ## Get coverage level and get protein annotations
    dict_edit, cov_dict, dict_prot_annot = get_edit_level_annotated(all_position_info_rna, result_file, unique_positions, rna_cov_threshold, edit_threshold)
    return dict_edit, cov_dict, dict_prot_annot


def call_RNA_editing_all(rna_dir, samples_list, ref_name, out_rna_edit, out_cov, out_rna_annot, out_unique_positions,
                         rna_cov_threshold=30, edit_threshold=0.03, n_cores=32):
    all_samples = samples_list

    unique_positions = set()

    ## get unique positions
    start_time = time.time()  # Start timing
    print(f"Computing RNA edit from {len(all_samples)} samples", flush=True)

    print("Searching for unique edit sites...", flush=True)
    if len(all_samples)>1:
        with Pool(n_cores) as pool:
            # First loop for editing
            all_dict_edit = list(tqdm(pool.starmap(compute_editing_sample, [(sample, rna_dir, ref_name, rna_cov_threshold, edit_threshold) for sample in all_samples]), total=len(all_samples)))

    else:
        all_dict_edit = [compute_editing_sample(all_samples[0], rna_dir, ref_name, rna_cov_threshold, edit_threshold)]
    # #### SEQUENTUAL. UNCOMMENT FOR DEBUGGING  ############################################
    # all_dict_edit = []
    # for sample in tqdm(all_samples):
    #     rna_edit_dict = compute_editing_sample(sample, rna_dir, ref_name, rna_cov_threshold, edit_threshold)
    #     all_dict_edit.append(rna_edit_dict)
    ########################################################################################
    for rna_edit_dict in all_dict_edit:
        if len(rna_edit_dict) > 0:
            unique_positions = unique_positions.union(set(rna_edit_dict.keys()))

    unique_positions = list(unique_positions)
    print(unique_positions)
    with open(out_unique_positions, 'w') as out:
        for pos in unique_positions:
            out.write(pos+'\n')
    print(f"Found {len(unique_positions)} sites.", flush=True)
    print(f"Execution time: {start_time - time.time()} seconds", flush=True)

    # # # Second loop for coverage
    # with Pool(n_cores) as pool:
    #     all_dict_edit, all_dict_cov, all_dict_prot_annot = list(tqdm(pool.starmap(compute_coverage_sample,
    #     [(sample, rna_dir, ref_name, rna_cov_threshold, edit_threshold, unique_positions) for sample in
    #         all_samples]), total=len(all_samples)))

    #### SEQUENTUAL. UNCOMMENT FOR DEBUGGING  ############################################
    start_time = time.time()
    print("Searching for coverage...", flush=True)
    all_dict_edit, all_dict_cov = [], []
    all_dict_prot_annot = dict()
    for sample in tqdm(all_samples):
        edit, cov, prot_annot = compute_coverage_sample(sample, rna_dir, ref_name, rna_cov_threshold, edit_threshold,
                                                        unique_positions)
        for pos in prot_annot.keys():
            if pos not in all_dict_prot_annot.keys():
                all_dict_prot_annot[pos] = prot_annot[pos]
        all_dict_edit.append(edit)
        all_dict_cov.append(cov)
    ########################################################################################
    # re-define unique positions
    unique_positions_updated = []
    for pos in unique_positions:
        if pos in all_dict_prot_annot.keys():
            unique_positions_updated.append(pos)
    print(f"{start_time - time.time()} seconds.", flush=True)

    ## Write unique positions to file


    ## Write to file
    start_time = time.time()
    print("Writing to file edits...", flush=True)
    write_edit_to_file(unique_positions_updated, all_dict_edit, out_rna_edit, all_samples)
    print(f"{start_time - time.time()} seconds.", flush=True)

    start_time = time.time()
    print("Writing to file coverage...", flush=True)
    write_cov_to_file(unique_positions_updated, all_dict_cov, out_cov, all_samples)
    print(f"{start_time - time.time()} seconds.", flush=True)

    start_time = time.time()
    print("Writing to file annotations...", flush=True)
    write_annot_to_file(all_dict_prot_annot, out_rna_annot)
    print(f"{start_time - time.time()} seconds.", flush=True)


# DNA_dir = './example/out_DNA/'
# RNA_dir = './example/out_RNA/'

# def annotate_gene(gff_file, pos):
#     return


def filter_final_RNA_editing(rna_edit_study, rna_cov_study, dna_edit_study, dna_cov_study, prot_annot_study, gff_file, out_dir, bacteria_name, cov_threshold):
    """
    """

    print(f"Starting filtering...")


    rna_df = pd.read_csv(rna_edit_study, index_col='position').T
    rna_df.columns = rna_df.columns.astype(str)

    print(f"Total potential RNA sites: {len(rna_df.columns)}")

    # potential_sites = [x for x in rna_df.columns if x not in all_dna_snps]
    # print(f"Number of filtered RNA sites: {len(potential_sites)}")

    # Filter out sites that had any editing in DNA
    dna_snp_df = pd.read_csv(dna_edit_study, index_col='position').T
    dna_snp_df.columns = dna_snp_df.columns.astype(str)

    filter1_list = []
    for site in rna_df.columns:
        if dna_snp_df[site].sum()==0:
            filter1_list.append(site)
    print(f"Number of sites passed DNA SNPs filter: {len(filter1_list)}")
    print(dna_snp_df)
    ## Filter out sites that are not present in DNA
    dna_cov_df = pd.read_csv(dna_cov_study, index_col='position').T
    print(dna_cov_study)
    dna_cov_df.columns = dna_cov_df.columns.astype(str)

    filtered_sites = []
    dna_cov = []
    dna_reads = []

    ### Read annotations for each site
    prot_anot_df = pd.read_csv(prot_annot_study, sep='\t')
    #print(prot_anot_df)
    anot_dict = {str(row['Position']): {col: row[col] for col in prot_anot_df.columns if col != 'Position'} for _, row in
                 prot_anot_df.iterrows()}
    for site in filter1_list:
        maxA = dna_cov_df[site+"_A"].max()
        maxC = dna_cov_df[site + "_C"].max()
        maxG = dna_cov_df[site + "_G"].max()
        maxT = dna_cov_df[site+"_T"].max()
        ref_base = anot_dict[site]['Old_base']
        #print(maxA, maxC, maxG, maxT)
        if ((ref_base=='A' and maxA>cov_threshold) or (ref_base=='T' and maxT>cov_threshold)) \
                        and maxC<cov_threshold \
                        and maxG<cov_threshold:

            N_A = (dna_cov_df[site + "_A"] != 0).sum()
            N_C = (dna_cov_df[site + "_C"] != 0).sum()
            N_G = (dna_cov_df[site + "_G"] != 0).sum()
            N_T = (dna_cov_df[site + "_T"] != 0).sum()

            reads_A = dna_cov_df[site + "_A"].sum()
            reads_C = dna_cov_df[site + "_C"].sum()
            reads_G = dna_cov_df[site + "_G"].sum()
            reads_T = dna_cov_df[site + "_T"].sum()


            ######### True edits - A to G; T to C
            if (ref_base == 'A' and N_C/N_A<0.1 and N_T/N_A<0.1 and N_G<3) or (ref_base == 'T' and N_G/N_T<0.1 and N_A/N_T<0.1 and N_C<3):
            # if the reference is A - allow for 10% C, and coverage in G < 2
            # if the reference is T - allow for 10% G, and coverage in C < 2

                filtered_sites.append(site)
                dna_cov.append((N_A, N_C, N_G, N_T))
                dna_reads.append((reads_A, reads_C, reads_G, reads_T))

    print(f"Number of sites that has coverage without SNPs and don't have coverage in C and G in DNA: {len(filtered_sites)}")
    if len(filtered_sites)==0:
        print("Didn't find any RNA edit")
        return 1
    ## Get the RNA coverage dictionary
    rna_cov_df = pd.read_csv(rna_cov_study, index_col='position').T
    rna_cov = []
    rna_reads = []
    for site in filtered_sites:
        N_A = (rna_cov_df[site + "_A"] != 0).sum()
        N_C = (rna_cov_df[site + "_C"] != 0).sum()
        N_G = (rna_cov_df[site + "_G"] != 0).sum()
        N_T = (rna_cov_df[site + "_T"] != 0).sum()

        reads_A = rna_cov_df[site + "_A"].sum()
        reads_C = rna_cov_df[site + "_C"].sum()
        reads_G = rna_cov_df[site + "_G"].sum()
        reads_T = rna_cov_df[site + "_T"].sum()


        rna_cov.append((N_A,N_C,N_G,N_T))
        rna_reads.append((reads_A, reads_C, reads_G, reads_T))

    filtered_df = rna_df[filtered_sites]
    # Sort by prevalence

    non_zero_counts = (filtered_df != 0).sum()

    # Sort columns by frequency in descending order
    sorted_columns = non_zero_counts.sort_values(ascending=False).index

    # Reindex DataFrame based on sorted columns
    filtered_df = filtered_df[sorted_columns]

    filtered_df.to_csv(f'{out_dir}/{bacteria_name}_filtered.csv', index=False)

    annotated_df = pd.DataFrame({"Position": filtered_sites,
                                 "N_edits_above_threshold": [(filtered_df[site] !=0).sum() for site in filtered_sites],
                                # "N_RNA_samples_A": [x[0] for x in rna_cov],
                                # "N_RNA_samples_C": [x[1] for x in rna_cov],
                                # "N_RNA_samples_G": [x[2] for x in rna_cov],
                                # "N_RNA_samples_T": [x[3] for x in rna_cov],
                                # "N_DNA_samples_A": [x[0] for x in dna_cov],
                                # "N_DNA_samples_C": [x[1] for x in dna_cov],
                                # "N_DNA_samples_G": [x[2] for x in dna_cov],
                                # "N_DNA_samples_T": [x[3] for x in dna_cov],
                                 "RNA_reads_A": [x[0] for x in rna_reads],
                                 "RNA_reads_C": [x[1] for x in rna_reads],
                                 "RNA_reads_G": [x[2] for x in rna_reads],
                                 "RNA_reads_T": [x[3] for x in rna_reads],
                                 "DNA_reads_A": [x[0] for x in dna_reads],
                                 "DNA_reads_C": [x[1] for x in dna_reads],
                                 "DNA_reads_G": [x[2] for x in dna_reads],
                                 "DNA_reads_T": [x[3] for x in dna_reads],
                                 })
    #pdb.set_trace()

    annotated_df = prot_anot_df.merge(annotated_df, on='Position', how='right')

    annotated_df['average_edit'] = annotated_df.apply(lambda row: get_edit_level_one(row['Old_base'],
                                                                                     row['RNA_reads_A'],
                                                                                     row['RNA_reads_C'],
                                                                                     row['RNA_reads_G'],
                                                                                     row['RNA_reads_T']), axis=1)


    annotated_df.to_csv(f'{out_dir}/{bacteria_name}_filtered_annotated.tsv', sep='\t', index=False)
    ## Annotate each site
    # with open(f'{out_dir}/{bacteria_name}_filtered_annotated.tsv', 'w') as out:
    #     out.write(f"ID\tRef\tRegion\tPosition\tN_RNA_samples_edit\tN_RNA_samples_A\tN_RNA_samples_C\tN_RNA_samples_G\tN_RNA_samples_T\tN_DNA_samples_A\tN_DNA_samples_C\tN_DNA_samples_G\tN_DNA_samples_T\tAvg_edit\tN_RNA_reads_A\tN_RNA_reads_C\tN_RNA_reads_G\tN_RNA_reads_T\tN_DNA_reads_A\tN_DNA_reads_C\tN_DNA_reads_G\tN_DNA_reads_T\t.\t.\t.\tGene_info\n")
    #     with open(gff_file, 'r') as f:
    #
    #         for line in f.readlines():
    #             line = line.strip('\n')
    #             if line[0]!='#': # skip comments
    #                 fields = line.split('\t')
    #                 start = int(fields[3])
    #                 end = int(fields[4])
    #                 for i, site in enumerate(filtered_sites):
    #                     if int(site)<end and int(site)>start:
    #                         nonzero = filtered_df[filtered_df[site]!=0][site]
    #                         freq = len(nonzero)
    #                         if freq>1:
    #                             avg_edit = nonzero.mean()
    #                             out.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{site}\t{freq}\t{rna_cov[i][0]}\t{rna_cov[i][1]}\t{rna_cov[i][2]}\t{rna_cov[i][3]}\t{dna_cov[i][0]}\t{dna_cov[i][1]}\t{dna_cov[i][2]}\t{dna_cov[i][3]}\t{avg_edit}\t{rna_reads[i][0]}\t{rna_reads[i][1]}\t{rna_reads[i][2]}\t{rna_reads[i][3]}\t{dna_reads[i][0]}\t{dna_reads[i][1]}\t{dna_reads[i][2]}\t{dna_reads[i][3]}\t{fields[5]}\t{fields[6]}\t{fields[7]}\t{fields[8]}\n")
    return 1


#################
## PARAMETERS
#################

# DNA_dir = './out_DNA/'
# RNA_dir = "/data/metaEdit/out_RNA/"
# REF_dir = '../../data/References/'
# ref_name = 'GCF_000008865.2_ASM886v2_genomic'
# bacteria_name = 'ecoli'
#
# samples_list = './data/samples.txt'
#
# rna_cov_threshold=30
# dna_cov_threshold=30
# edit_threshold=0.03
#
# out_rna_dir = f'./out_rna_edit_cov{dna_cov_threshold}_1015/'
# out_dna_dir = f'./out_dna_snps_cov{dna_cov_threshold}_1015/'

#################

# os.makedirs(out_rna_dir, exist_ok=True)
# os.makedirs(out_dna_dir, exist_ok=True)
#
# gff_file = REF_dir + f'{bacteria_name}/{ref_name}.gff'
#
# rna_edit_study = f'./{out_rna_dir}/{bacteria_name}_edit.csv'
# rna_cov_study = f'./{out_rna_dir}/{bacteria_name}_cov.csv'
#
# dna_edit_study = f'./{out_dna_dir}/{bacteria_name}_snps.csv'
# dna_cov_study = f'./{out_dna_dir}/{bacteria_name}_cov.csv'
#
# print("Helo")
#
# call_RNA_editing_all(RNA_dir, samples_list, ref_name, rna_edit_study,
#                                                                 rna_cov_study,
#                                                                 rna_cov_threshold,
#                                                                 edit_threshold)
#
# get_dna_snps_whole_study(DNA_dir, samples_list, ref_name, rna_edit_study, dna_edit_study,
#                          dna_cov_study, cov_threshold=dna_cov_threshold)
#
#
# =filter_final_RNA_editing(rna_edit_study, dna_edit_study, dna_cov_study, gff_file, out_rna_dir, bacteria_name, cov_threshold=dna_cov_threshold)

def main(DNA_dir, RNA_dir, REF_dir, ref_name, bacteria_name,
         rna_cov_threshold, dna_cov_threshold, edit_threshold,
         out_rna_dir, out_dna_dir, n_cores=32, rna_only=False, DNA_sample=None, samples_list=None, sample=None):
    print("Starting filtering...")
    os.makedirs(out_rna_dir, exist_ok=True)
    os.makedirs(out_dna_dir, exist_ok=True)

    if bacteria_name==DNA_sample:
        bacteria_name_prefix = ""
    else:
        bacteria_name_prefix = bacteria_name

    gff_file = REF_dir + f'{bacteria_name_prefix}/{ref_name}.gff'
    if not os.path.exists(gff_file):
        gff_file = REF_dir + f'{bacteria_name_prefix}/{ref_name}.gff3'

    rna_edit_study = f'{out_rna_dir}/{bacteria_name}_edit.csv'
    rna_cov_study = f'{out_rna_dir}/{bacteria_name}_cov.csv'
    rna_annot_study = f'{out_rna_dir}/{bacteria_name}_annot.tsv'
    out_unique_positions = f'{out_rna_dir}/{bacteria_name}_unique_pos.txt'

    if sample is not None:
        samples_list = [sample]
    elif samples_list is not None:
        samples_list = [x.strip('\n') for x in open(samples_list).readlines()]

    if DNA_sample is not None:
        dna_edit_study = f'{out_dna_dir}/{bacteria_name}_snps-{DNA_sample}.csv'
        dna_cov_study = f'{out_dna_dir}/{bacteria_name}_cov-{DNA_sample}.csv'
        samples_list = [DNA_sample]
    else:
        dna_edit_study = f'{out_dna_dir}/{bacteria_name}_snps.csv'
        dna_cov_study = f'{out_dna_dir}/{bacteria_name}_cov.csv'
        call_RNA_editing_all(RNA_dir, samples_list, ref_name, rna_edit_study,
                             rna_cov_study, rna_annot_study, out_unique_positions,
                             rna_cov_threshold, edit_threshold, n_cores=n_cores)

    if not rna_only:
        get_dna_snps_whole_study(DNA_dir, samples_list, ref_name, rna_edit_study, dna_edit_study,
                                 dna_cov_study, cov_threshold=dna_cov_threshold, n_cores=n_cores, out_dir = out_dna_dir)

        filter_final_RNA_editing(rna_edit_study, rna_cov_study, dna_edit_study, dna_cov_study, rna_annot_study, gff_file, out_rna_dir, bacteria_name, cov_threshold=dna_cov_threshold)


if __name__ == '__main__':
    start_time = time.time()  # Start timing

    parser = argparse.ArgumentParser(description="Command line arguments for the script")

    parser.add_argument("DNA_dir")
    parser.add_argument("RNA_dir")
    parser.add_argument("REF_dir")
    parser.add_argument("ref_name")
    parser.add_argument("bacteria_name")
    parser.add_argument("rna_cov_threshold", type=int)
    parser.add_argument("dna_cov_threshold", type=int)
    parser.add_argument("edit_threshold", type=float)
    parser.add_argument("out_rna_dir")
    parser.add_argument("out_dna_dir")
    parser.add_argument("cores", type=int)
    parser.add_argument("--samples_list")
    parser.add_argument("--rna_only", action="store_true", default=False)
    parser.add_argument("--DNA_sample")
    parser.add_argument("--sample")

    args = parser.parse_args()

    main(args.DNA_dir, args.RNA_dir, args.REF_dir, args.ref_name,
         args.bacteria_name, args.rna_cov_threshold,
         args.dna_cov_threshold, args.edit_threshold, args.out_rna_dir,
         args.out_dna_dir, args.cores, args.rna_only, sample=args.sample, DNA_sample=args.DNA_sample,  samples_list=args.samples_list)

    print(f"Execution time: {time.time() - start_time} seconds")
    # DNA_dir = './example/out_DNA/'
    # RNA_dir = "./example/out_RNA/"
    # REF_dir = '../../data/References/'


    # ref_name = 'GCF_000008865.2_ASM886v2_genomic'
    # bacteria_name = 'ecoli'

    # ref_name = 'GCF_000154465.1_ASM15446v1_genomic'
    # bacteria_name = 'Alistipes_putredinis'
    #
    # ref_name = 'GCF_008572305.1_ASM857230v1_genomic'
    # bacteria_name = 'Bacteroides_caccae'
    #
    # ref_name = 'GCF_001406095.1_13414_6_48_genomic'
    # bacteria_name = 'Bacteroides_uniformis'
    #
    # ref_name = 'GCF_000012825.1_ASM1282v1_genomic'
    # bacteria_name = 'Bacteroides_vulgatus'
    #
    # ref_name = 'Dialister_invisus'
    # bacteria_name = 'Dialister_invisus'
    #
    # ref_name = 'GCF_003312465.1_ASM331246v1_genomic'
    # bacteria_name = 'Faecalibacterium_prausnitzii'
    #
    # # ref_name = 'GCF_000307495.1_Para_merd_CL09T00C40_V1_genomic'
    # # bacteria_name = 'Parabacteroides_merdae'
    #
    # samples_list = './data/samples.txt'
    #
    # rna_cov_threshold=30
    # dna_cov_threshold=30
    # edit_threshold=0.03
    #
    # # out_rna_dir = f'./out_rna_edit_cov{dna_cov_threshold}_1015_tmp/'
    # # out_dna_dir = f'./out_dna_snps_cov{dna_cov_threshold}_1015_tmp/'
    #
    # out_rna_dir = f'./out_RNA_all_1016/'
    # out_dna_dir = f'./out_DNA_all_1016/'
    #
    # main(DNA_dir, RNA_dir, REF_dir, ref_name,
    #      bacteria_name, samples_list, rna_cov_threshold,
    #      dna_cov_threshold, edit_threshold, out_rna_dir,
    #      out_dna_dir)
    #
