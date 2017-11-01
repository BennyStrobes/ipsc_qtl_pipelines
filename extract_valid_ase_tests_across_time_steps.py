import numpy as np
import os
import sys
import pdb

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


# Limit to rows (sites) that pass our allelic counts filters
def determine_allelic_count_rows_that_pass_filters(filtered_counts, valid_columns, min_reads, min_samples_per_time_step, min_fraction_biallelic):
    valid_indices = []  # Keep track of rows that pass allelic filters
    num_rows, num_cols = filtered_counts.shape  # Get dimensions of matrix
    # Loop through all of rows (all of heterozygous sites)
    for row_num in range(num_rows):
        pass_filter = True  # Initialize pass_filter to True
        # Must pass filter for samples from ALL of time steps
        for time_step in range(16):
            time_step_columns = valid_columns[time_step]
            time_step_allelic_counts = np.squeeze(np.asarray(filtered_counts[row_num,time_step_columns]))
            # First check if there are at least min_samples_per_time_step samples with heterozygous sites
            if np.sum(time_step_allelic_counts != 'Nan') < min_samples_per_time_step:
                pass_filter = False
                continue
            # Compute allelic imbalence for each of the sampels
            allelic_imbalences = []
            for count_string in time_step_allelic_counts:
                if count_string == 'Nan': # Ignore nan values
                    continue
                ref_count = float(count_string.split('_')[0])  # extract reference allele counts
                total_count = float(count_string.split('_')[1])  # Extract total counts
                if total_count >= min_reads:
                    allelic_imbalences.append(abs((ref_count/total_count) - .5))
                # If there are fewer than min_reads total counts, assume sample is not biallelic
                else:
                    allelic_imbalences.append(.5)
            # now check if at least min_fraction_biallelic of the heterozygous samples show bi-allelic expression
            fraction_biallelic = np.sum(np.asarray(allelic_imbalences) < .49)/len(allelic_imbalences)
            if fraction_biallelic < min_fraction_biallelic:
                pass_filter = False
                continue
        # Check to see if row has passed filter for all 16 time steps
        if pass_filter == True:
            valid_indices.append(row_num)
    return valid_indices


# Extract vector of length number of time steps where each element of the vector is a matrix describing the total counts matrix (sitesXnum_samples) for that time step
# Extract vector of length number of time steps where each element of the vector is a matrix describing the allelic imbalence  matrix (sitesXnum_samples) for that time step
def extract_allelic_count_information(allelic_counts_file, chrom_num, min_reads, min_samples_per_time_step, min_fraction_biallelic, time_step):
    ############ Load in data line by line
    f = open(allelic_counts_file)
    head_count = 0  # To skip header
    site_ids = []  # Initialize array to keep ordered list of heterozygous site ids
    positions = []  # Initialize array to keep ordered list of heterozygous site positions
    all_counts = []  # Initialize array to keep ordered list of allelic_count arrays at each site
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Header
            head_count = head_count + 1
            all_samples = data[1:]  # Ordered array of sample names
            continue
        # Extract relevent data
        site_id = data[0]
        counts = data[1:]
        # Limit to sites on this chromosome
        site_chrom = int(site_id.split('_')[0])
        if site_chrom != chrom_num:  # skip sites not on this chromosome
            continue
        # Keep track of sites on this chromsome
        site_ids.append(site_id)
        positions.append(int(site_id.split('_')[1]))
        all_counts.append(counts)
    # put in correct format
    filtered_counts = np.asmatrix(all_counts)  # convert from list of arrays to a matrix
    all_samples = np.asarray(all_samples)
    site_ids = np.asarray(site_ids)

    ############# Get indices of samples for each time step
    valid_columns = []  # Initialize list of length number of time steps, where each element is an ordered array that contains the sample indices corresponding to that time step
    for temp_time_step in range(16):  # Loop through 16 time steps
        indices = []  # keep track of indices of samples in this time step
        for index, sample_id in enumerate(all_samples):
            sample_info = sample_id.split('_')
            cell_line = sample_info[0]
            sample_time_step = int(sample_info[1])
            if sample_time_step != temp_time_step:
                continue
            # If we've made it this far then the sample is at the correct time step
            indices.append(index)
        valid_columns.append(np.asarray(indices))

    # Limit to rows (sites) that pass our allelic counts filters
    rows_that_pass_allelic_filters = determine_allelic_count_rows_that_pass_filters(filtered_counts, valid_columns, min_reads, min_samples_per_time_step, min_fraction_biallelic)
    # Now apply those filters
    filtered_counts = filtered_counts[rows_that_pass_allelic_filters, :]
    positions = np.asarray(positions)[rows_that_pass_allelic_filters]
    site_ids = np.asarray(site_ids)[rows_that_pass_allelic_filters]

    # Filter to correct columns (samples at this time step)
    correct_columns = valid_columns[time_step]
    filtered_counts = filtered_counts[:,correct_columns]


    # Convert from matrix of refCounts_totalCounts (strings) to matrix of allelic imbalence
    allelic_imbalence_matrix, total_counts_mat = create_allelic_imbalence_matrix(filtered_counts)
    return allelic_imbalence_matrix, total_counts_mat, positions, site_ids


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
        #site_id = data[0] + '_' + data[1] + '_' + data[3] + '_' + data[4]
        site_id = data[2]
        # Append matrices
        reg_prob_matrix.append(line_probs)
        reg_site_positions.append(position)
        reg_site_ids.append(site_id)
    reg_prob_matrix = np.asmatrix(reg_prob_matrix)
    return reg_prob_matrix, np.asarray(reg_site_positions), np.asarray(reg_site_ids)

# Extract list of cell lines that are found in all time steps (some time steps have different numbers of cell lines)
# When calling maf cutoff, we are to use only cell lines found in all time steps
def extract_list_of_cell_lines_in_all_time_steps(input_file):
    # First, extract array of sample ids from quantile normalized expression matrix
    f = open(input_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            sample_ids = data[1:]
            continue
        continue
    f.close()
    # Now loop through all time steps to determine which cell lines were observed in which time step
    cell_lines_all_time_steps = []
    for time_step in range(0,16):
        temp_dicti = []
        for sample_id in sample_ids:
            # check if sample_id is from the current iteration's time step
            if sample_id.endswith('_' + str(time_step)) == False:
                continue
            # Get cell line from sample id
            cell_line = sample_id.split('_')[0]
            temp_dicti.append(cell_line)
        cell_lines_all_time_steps.append(temp_dicti)
    return cell_lines_all_time_steps

#  Extract of array of length number of time steps, where each element of array is another array containing indices of samples
def get_valid_samples_across_time_steps(total_counts_matrices, het_site_index, min_reads):
    valid_samples = []
    #  Loop throug htime steps
    for time_step in range(16):
        # Get total counts matrix for this time step
        time_step_total_counts = total_counts_matrices[time_step]
        # Limit to specific heterozygous site of interest
        total_counts_vector = time_step_total_counts[het_site_index,:]
        # Get samples that have reads mapping
        valid_samples.append(total_counts_vector >= min_reads)
    return valid_samples

# Find indices of reg_site_positions that are within distance of het_site_position
def find_reg_site_positions_within_distance_cis(het_site_position, reg_site_positions, distance):
    # Compute absolute distance from each reg_site and the current heterozygous site
    absolute_difference_vector = abs(reg_site_positions - het_site_position)
    # Find all indices where distance between regulatory site and heterozygous site is less than distance
    valid_indices = np.where(absolute_difference_vector <= distance)[0]
    return valid_indices

def get_heterozygous_maf(reg_site_vector):
    num_heterozygous = np.sum(reg_site_vector)
    num_homozygous = len(reg_site_vector) - num_heterozygous
    total_samples = num_heterozygous + num_homozygous
    return min(num_heterozygous/total_samples, num_homozygous/total_samples)

# Make sure maf >= $min_fraction_in_test_group for EVERY time step
def pass_maf_filter_all_time_steps(reg_prob_matrices, reg_site_index, valid_samples_across_time_steps,min_fraction_in_test_group):
    pass_filter = True  # Initialize to be true

    # Now check if false in any of time steps
    for temp_time_step in range(16):
        reg_prob_mat = reg_prob_matrices[temp_time_step]  # get reg. probability matrix for current time step
        valid_samples = valid_samples_across_time_steps[temp_time_step]  # Get valid samples (samples that have biallelic expression at heterozygous site) for this time step
        reg_site_vector = np.squeeze(np.asarray(reg_prob_mat[reg_site_index,:][:,valid_samples]))  # Get vector form
        af = get_heterozygous_maf(reg_site_vector)  # compute alllele frequence
        if af < min_fraction_in_test_group:
            pass_filter = False
    return pass_filter


# Main driver function
def driver(het_prob_genotype_file, allelic_counts_file, output_file, min_reads, min_samples_per_time_step, min_fraction_in_test_group, distance, chrom_num, min_fraction_biallelic):
    #####################################
    # Preprocess the data
    #####################################
    # Extract vector of length number of time steps where each element of the vector is a matrix describing the total counts matrix (sitesXnum_samples) for that time step
    # Extract vector of length number of time steps where each element of the vector is a matrix describing the allelic imbalence  matrix (sitesXnum_samples) for that time step
    total_counts_matrices = []
    allelic_imbalence_matrices = []
    for temp_time_step in range(16):
        print(temp_time_step)
        # For one time step
        allelic_imbalence_mat, total_counts_mat, het_site_positions, het_site_ids = extract_allelic_count_information(allelic_counts_file, chrom_num, min_reads, min_samples_per_time_step, min_fraction_biallelic, temp_time_step)
        total_counts_matrices.append(total_counts_mat)
        allelic_imbalence_matrices.append(allelic_imbalence_mat)

    # Extract list of length number of time steps where each element is a dictionary that contains the cell lines observed for that time step
    # When calling maf cutoff, we require it to pass the maf in each observed time step
    cell_lines_all_time_steps = extract_list_of_cell_lines_in_all_time_steps(allelic_counts_file)

    # Extract heterozygous probabilities of all regulatory variants on this chromosome (do this for each of 16 time steps)
    # Also, orders columns in correct order with respect to samples in allelic_imbalence_matrices
    reg_prob_matrices = []
    for temp_time_step in range(16):
        print(temp_time_step)
        reg_prob_matrix, reg_site_positions, reg_site_ids = extract_regulatory_variant_information(het_prob_genotype_file, cell_lines_all_time_steps[temp_time_step], chrom_num)
        reg_prob_matrices.append(reg_prob_matrix)

    # Initialize Output handle
    t = open(output_file, 'w')

    num_het_sites = len(het_site_positions)

    # Loop through all instances of the dependent variable (allelic imbalence)
    for het_site_index in range(num_het_sites):
        # Extract quantities of allelic_imbalence related to this iterations het_site_index
        het_site_position = het_site_positions[het_site_index]
        het_site_id = het_site_ids[het_site_index]
        het_site_rs_id = het_site_id.split('_')[2]

        #  Extract of array of length number of time steps, where each element of array is another array containing indices of samples
        valid_samples_across_time_steps = get_valid_samples_across_time_steps(total_counts_matrices, het_site_index, min_reads)

        # Find indices of reg_site_positions that are within distance of het_site_position
        reg_site_indices = find_reg_site_positions_within_distance_cis(het_site_position, reg_site_positions, distance)
        
        # Now loop through all valid regulatory sites
        for reg_site_index in reg_site_indices:
            # Extract quantities of the regulatory site related to this iterations reg_site_index
            reg_site_position = reg_site_positions[reg_site_index]
            reg_site_id = reg_site_ids[reg_site_index]

            # Make sure maf >= $min_fraction_in_test_group for EVERY time step
            pass_filter = pass_maf_filter_all_time_steps(reg_prob_matrices, reg_site_index, valid_samples_across_time_steps,min_fraction_in_test_group)
            if pass_filter == False:
                continue
            # Make sure regulatory variant is not the heterozygous variant!
            if reg_site_id == het_site_rs_id:
                continue

            # Limit to variants for which we have a known rsid
            if reg_site_id == '.' or het_site_rs_id == '.':
                continue

            # We have now passed all filters
            # Ie, this variant gene pair is a valid test
            t.write(het_site_id + '\t' + reg_site_id + '\n')
    t.close()







###############################################################################
# Objective
############
# The goal of this script is to create a list of variant-gene pairs that pass all
# of our ase-qtl filters.
# The filters are as follows:
# A variant gene pair must pass all of the following in EVERY time step:
#### 1. distance between reg. variant and heterozygous site is < distance
#### 2. There are at least min_samples_per_time_step that are heterozygous at the exonic variant
#### 3. Of the the heterozygous samples, at least min_fraction_biallelic % of them show biallelic expression
####- where biallelic expression is both: a) greater than min_reads total reads and b) allelic imbalence < .49
#### 4. Regulatory variant must have MAF >= min_fraction_in_test_group
###############################################################################



##########################################
# Input Data
##########################################
het_prob_genotype_file = sys.argv[1]
allelic_counts_file = sys.argv[2]
output_file = sys.argv[3]
min_reads = int(sys.argv[4])  # Minimum number of reads on heterozygous site to consider a sample valid
min_samples_per_time_step = int(sys.argv[5])  # Minimum number of valid samples we will run a qtl test on
min_fraction_in_test_group = float(sys.argv[6])  # Minimum fraction of valid samples with less popular version of regulatory variant (homozygous reg variant vs heterozygous reg variant)
distance = float(sys.argv[7])  # Cis eqtl distance (EAGLE used 200 KB)
chrom_num = int(sys.argv[8])
min_fraction_biallelic = float(sys.argv[9])  # Minimum fraction of heterozygous samples that show biallelic expression


driver(het_prob_genotype_file, allelic_counts_file, output_file, min_reads, min_samples_per_time_step, min_fraction_in_test_group, distance, chrom_num, min_fraction_biallelic)
