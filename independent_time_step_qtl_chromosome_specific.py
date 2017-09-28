import numpy as np
import os
import sys
import pdb
from ase_qtl import run_ase_qtl

##########################################
# Functions
##########################################
# Convert from matrix of refCounts_totalCounts (strings) to matrix of allelic imbalence
def create_allelic_imbalence_matrix(filtered_counts):
    row_num, col_num = filtered_counts.shape
    # Initialize allelic imbalence mat
    allelic_imbalence_mat = np.zeros((row_num,col_num))
    total_counts_mat = np.zeros((row_num,col_num))
    # loop through rows and columns
    for ii in range(row_num):
        for jj in range(col_num):
            if filtered_counts[ii,jj] == 'Nan':  # Cell line, site is not heterozygous
                allelic_imbalence_mat[ii,jj] = np.nan
                total_counts_mat[ii,jj] = np.nan
            else: # Cell line, site is heterozygous
                # Extract ref counts and total counts from string
                info = filtered_counts[ii,jj].split('_')
                ref_count = float(info[0])
                total_count = float(info[1])
                total_counts_mat[ii,jj] = total_count
                # compute allelic imbalence
                if total_count > 0:
                    allelic_imbalence = abs((ref_count/total_count) - .5)
                    allelic_imbalence_mat[ii,jj] = allelic_imbalence
                else:
                    allelic_imbalence_mat[ii,jj] = np.nan
    return allelic_imbalence_mat, total_counts_mat

# Extract matrix of allelic imbalence for all sites on this chromosome
def extract_allelic_count_information(allelic_counts_file, chrom_num, time_step):
    ############ Load in data
    data = np.loadtxt(allelic_counts_file, dtype=str)
    # All possible heterozygous sites (across chromosomes)
    all_sites = data[1:,0]
    # All possible samples (across time steps)
    all_samples = data[0,1:]
    # all count data:
    all_counts = data[1:,1:]
    ############# Filter to sites on this chromosome
    positions = []  # Initialize array to keep track of positions of heterozygous sites on this chromosome
    site_ids = [] # Initialize array to keep track of names of heterozygous sites being used
    valid_rows = []  # Keep track of row indices of all_sites that are on chrom_num
    # Loop through all sites
    for index, site in enumerate(all_sites):
        info = site.split('_')
        site_chrom = int(info[0])
        if site_chrom != chrom_num:
            continue
        # If we've made it this far then site is on the correct chromosome
        position = int(info[1])
        valid_rows.append(index)
        positions.append(position)
        site_ids.append(site)
    ############# Filter to samples at the correct time step
    ordered_cell_lines = [] # Initialize array to keep track of cell lines
    valid_columns = []  # Keep track of column indices of valid samples
    for index, sample_id in enumerate(all_samples):
        sample_info = sample_id.split('_')
        cell_line = sample_info[0]
        sample_time_step = int(sample_info[1])
        if sample_time_step != time_step:
            continue
        # If we've made it this far then the sample is at the correct time step
        ordered_cell_lines.append(cell_line)
        valid_columns.append(index)
    # Filter to correct rows (only sites on this chromosome) and correct columns (samples at this time step)
    filtered_counts = all_counts[valid_rows,:][:,valid_columns]
    # Convert from matrix of refCounts_totalCounts (strings) to matrix of allelic imbalence
    allelic_imbalence_matrix, total_counts_mat = create_allelic_imbalence_matrix(filtered_counts)
    return allelic_imbalence_matrix, total_counts_mat, positions, site_ids, ordered_cell_lines

# Extract heterozygous probabilities of all regulatory variants on this chromosome
def extract_regulatory_variant_information(het_prob_genotype_file, ordered_cell_lines, chrom_num):
    reg_prob_matrix = []  # Initialize array to keep track of reg probs
    reg_site_positions = []  # Initialize array to keep track of all positions of regulatory sites on this chromosome
    reg_site_ids = []  # Initialize array to keep track of names of regulatory sites on this chromosome
    # loop through het_prob_genotype_file
    f = open(het_prob_genotype_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if line.startswith('#CHROM'):  # header with sample ids
            sample_indices = []  # ordered list of indices corresponding to cell lines we have rna-seq for at this time step
            # Make sure genotype and ae data is in same order
            for i,cell_line in enumerate(ordered_cell_lines):
                for index, val in enumerate(data):
                    if val == cell_line:
                        sample_indices.append(index)
            if np.array_equal(np.asarray(data)[sample_indices], np.asarray(ordered_cell_lines)) == False: # Simple error checking
                print('FATAL ERROR: GENOTYPE Cell lines doesnt match ordered_cell_lines')
                pdb.set_trace()
            continue
        if line.startswith('#'):  # ignore information headers of vcf file
            continue
        # Normal data
        line_chrom_num = int(data[0])
        if line_chrom_num != chrom_num:  # Ignore lines not on the correct chromosome
            continue
        # Het. probs for all lines we have rna seq for
        line_probs = np.asarray(data)[sample_indices].astype(float)
        # Position of variant
        position = int(data[1])
        # Name of variant
        site_id = data[0] + '_' + data[1] + '_' + data[3] + '_' + data[4]
        # Append matrices
        reg_prob_matrix.append(line_probs)
        reg_site_positions.append(position)
        reg_site_ids.append(site_id)
    reg_prob_matrix = np.asmatrix(reg_prob_matrix)
    return reg_prob_matrix, reg_site_positions, reg_site_ids

# Main driver function
def driver(het_prob_genotype_file, allelic_counts_file, output_file_prefix, min_reads, min_samples_per_time_step, min_fraction_in_test_group, distance, statistical_test, normalization_method, chrom_num, time_step, min_fraction_biallelic):
    #####################################
    # Preprocess the data
    #####################################
    # Extract matrix of allelic imbalence for all sites on this chromosome and all samples in this time step
    allelic_imbalence_mat, total_counts_mat, het_site_positions, het_site_ids, ordered_cell_lines = extract_allelic_count_information(allelic_counts_file, chrom_num, time_step)
    
    # Extract heterozygous probabilities of all regulatory variants on this chromosome
    reg_prob_matrix, reg_site_positions, reg_site_ids = extract_regulatory_variant_information(het_prob_genotype_file, ordered_cell_lines, chrom_num)
    
    #####################################
    # Run independent_time_step_qtl.py
    #####################################
    # Method to run aseqtl analysis on one chromosome
    run_ase_qtl(output_file_prefix, allelic_imbalence_mat, total_counts_mat, het_site_positions, het_site_ids, reg_prob_matrix, reg_site_positions, reg_site_ids, min_reads, min_samples_per_time_step, min_fraction_in_test_group, distance, statistical_test, normalization_method, chrom_num, min_fraction_biallelic)






##########################################
# Input Data
##########################################

het_prob_genotype_file = sys.argv[1]
allelic_counts_file = sys.argv[2]
output_file_prefix = sys.argv[3]
min_reads = int(sys.argv[4])  # Minimum number of reads on heterozygous site to consider a sample valid
min_samples_per_time_step = int(sys.argv[5])  # Minimum number of valid samples we will run a qtl test on
min_fraction_in_test_group = float(sys.argv[6])  # Minimum fraction of valid samples with less popular version of regulatory variant (homozygous reg variant vs heterozygous reg variant)
distance = float(sys.argv[7])  # Cis eqtl distance (EAGLE used 200 KB)
statistical_test = sys.argv[8]  # Statistical test to use in evaluating aseQTL (options are "wilcoxon" and "linear_regression")
normalization_method = sys.argv[9]  # Normalization method (only applicable if statistical_test="linear_regression", otherwise will be ignored). Options available are "standardize" and "quantile_normalize"
chrom_num = int(sys.argv[10])
time_step = int(sys.argv[11])
min_fraction_biallelic = float(sys.argv[12])  # Minimum fraction of heterozygous samples that show biallelic expression

# Main driver function
driver(het_prob_genotype_file, allelic_counts_file, output_file_prefix, min_reads, min_samples_per_time_step, min_fraction_in_test_group, distance, statistical_test, normalization_method, chrom_num, time_step, min_fraction_biallelic)

