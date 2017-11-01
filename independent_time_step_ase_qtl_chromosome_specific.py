import numpy as np
import os
import sys
import pdb
from scipy.stats import ranksums
from scipy.stats import rankdata
from scipy import stats
import statsmodels.api as sm


# We previously computed a list of variant-gene pairs that passed required filters
# Now extract dictionaries containing that information
def extract_valid_tests(valid_tests_input_file):
    # Initialize dictionaries
    valid_het_sites = {}
    valid_reg_sites = {}
    valid_tests = {}
    f = open(valid_tests_input_file)
    # Stream lines of input file
    for line in f:
        line = line.rstrip()
        data = line.split()
        # Extract relevent information
        het_site = data[0]
        reg_site = data[1]
        test_id = het_site + '_' + reg_site
        # Add relevent information to dictionaries
        valid_het_sites[het_site] = 1
        valid_reg_sites[reg_site] = 1
        valid_tests[test_id] = 1
    return valid_het_sites, valid_reg_sites, valid_tests



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
    return np.asmatrix(allelic_imbalence_mat), np.asmatrix(total_counts_mat)



def extract_allelic_count_information(allelic_counts_file, chrom_num, time_step, min_reads, valid_het_sites):
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
        # Limit to sites that are considered valid
        if site_id not in valid_het_sites:
            continue
        counts = data[1:]
        # Keep track of sites on this chromsome
        site_ids.append(site_id)
        positions.append(int(site_id.split('_')[1]))
        all_counts.append(counts)
    # put in correct format
    filtered_counts = np.asmatrix(all_counts)  # convert from list of arrays to a matrix
    all_samples = np.asarray(all_samples)
    site_ids = np.asarray(site_ids)

    ############# Get indices of samples for this time step
    correct_columns = []  # keep track of indices of samples in this time step
    for index, sample_id in enumerate(all_samples):
        sample_info = sample_id.split('_')
        cell_line = sample_info[0]
        sample_time_step = int(sample_info[1])
        if sample_time_step != time_step:
            continue
        # If we've made it this far then the sample is at the correct time step
        correct_columns.append(index)


    # Filter to correct columns (samples at this time step)
    filtered_counts = filtered_counts[:,correct_columns]


    # Convert from matrix of refCounts_totalCounts (strings) to matrix of allelic imbalence
    allelic_imbalence_matrix, total_counts_mat = create_allelic_imbalence_matrix(filtered_counts)
    return allelic_imbalence_matrix, total_counts_mat, positions, site_ids

# Extract list of all cell lines used in this time step
def extract_list_of_cell_lines_in_all_time_steps(input_file, time_step):
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
    cell_lines = []
    for sample_id in sample_ids:
        # check if sample_id is from the current iteration's time step
        if sample_id.endswith('_' + str(time_step)) == False:
            continue
        # Get cell line from sample id
        cell_line = sample_id.split('_')[0]
        cell_lines.append(cell_line)
    return cell_lines


# Extract heterozygous probabilities of all regulatory variants on this chromosome
def extract_regulatory_variant_information(het_prob_genotype_file, ordered_cell_lines, chrom_num, valid_reg_sites):
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

        # Only consider sites considered to be valid
        site_id = data[2]
        if site_id not in valid_reg_sites:
            continue
        # Append matrices
        reg_prob_matrix.append(line_probs)
        reg_site_positions.append(position)
        reg_site_ids.append(site_id)
    reg_prob_matrix = np.asmatrix(reg_prob_matrix)
    return reg_prob_matrix, np.asarray(reg_site_positions), np.asarray(reg_site_ids)

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

# Run ase-qtl for one regulatory variatnt - heterozygous site pair
def run_statistical_test(reg_site_vector, allelic_imbalence_vector_filtered_samples, statistical_test, normalization_method):
    # Simple error checking
    if len(reg_site_vector) != len(allelic_imbalence_vector_filtered_samples):
        print('Fatal error in ase-qtl test')
        pdb.set_trace()
    # Evaluate significance with wilcoxon-rank sum test
    if statistical_test == 'wilcoxon':
        # Get indices of all samples that are heterozygous at the regulatory locus
        het_indices = np.where(reg_site_vector >= .5)[0]
        # Get indices of all samples that are homozygous at the regulatory locus
        homo_indices = np.where(reg_site_vector < .5)[0]
        # Now create two vectors of allelic imbalence from the two groups
        het_ae = allelic_imbalence_vector_filtered_samples[het_indices]
        homo_ae = allelic_imbalence_vector_filtered_samples[homo_indices]
        statistic, pvalue = ranksums(het_ae, homo_ae)
    # Evaluate significance using linear regression
    elif statistical_test == 'linear_regression':
        if normalization_method == 'standardize':
            # Just standardize ae data
            normalized_ae = (allelic_imbalence_vector_filtered_samples - np.mean(allelic_imbalence_vector_filtered_samples))/np.std(allelic_imbalence_vector_filtered_samples)
        elif normalization_method == 'quantile_normalize':
            # Project ae data onto quantiles of a gaussian
            ranks = rankdata(allelic_imbalence_vector_filtered_samples)/(len(allelic_imbalence_vector_filtered_samples) + 1)
            normalized_ae = stats.norm.ppf(ranks)
        # Run linear regression
        reg_site_vector = sm.add_constant(reg_site_vector)  # Fit intercept
        fit = sm.OLS(normalized_ae, reg_site_vector).fit()
        pvalue = fit.pvalues[1]
        statistic = fit.params[1]
    return statistic, pvalue

# Main driver function
def driver(het_prob_genotype_file, allelic_counts_file, output_file, min_reads, distance, statistical_test, normalization_method, chrom_num, time_step, valid_het_sites, valid_reg_sites, valid_tests):
    #####################################
    # Preprocess the data
    #####################################
    # Extract vector of length number of time steps where each element of the vector is a matrix describing the total counts matrix (sitesXnum_samples) for that time step
    # Extract vector of length number of time steps where each element of the vector is a matrix describing the allelic imbalence  matrix (sitesXnum_samples) for that time step
    allelic_imbalence_mat, total_counts_mat, het_site_positions, het_site_ids = extract_allelic_count_information(allelic_counts_file, chrom_num, time_step, min_reads, valid_het_sites)

    # Extract list of all cell lines used in this time step
    cell_lines = extract_list_of_cell_lines_in_all_time_steps(allelic_counts_file, time_step)

    # Extract heterozygous probabilities of all regulatory variants on this chromosome for this time step
    # Also, orders columns in correct order with respect to samples in allelic_imbalence_matrices
    reg_prob_mat, reg_site_positions, reg_site_ids = extract_regulatory_variant_information(het_prob_genotype_file, cell_lines, chrom_num, valid_reg_sites)


    # Initialize output handle
    t = open(output_file, 'w')
    # Write header
    t.write('chom_num\thet_site_id\thet_site_position\tregulatory_site_id\tregulatory_site_position\tregulatory_allele_frequency\ttest_statistic\tpvalue\n')


    num_het_sites, num_samples = allelic_imbalence_mat.shape

    # First loop through all instances of the dependent variable (allelic imbalence)
    for het_site_index in range(num_het_sites):
        
        # Extract quantities of allelic_imbalence related to this iterations het_site_index
        allelic_imbalence_vector = np.squeeze(np.asarray(allelic_imbalence_mat[het_site_index,:]))
        total_counts_vector = np.squeeze(np.asarray(total_counts_mat[het_site_index,:]))
        het_site_position = het_site_positions[het_site_index]
        het_site_id = het_site_ids[het_site_index]
        het_site_rs_id = het_site_id.split('_')[2]

        # consider all samples that have het. site and have at least min_reads read mapping to the het site
        valid_samples = total_counts_vector >= min_reads

        # Filter allelic count vectors to valid_samples
        allelic_imbalence_vector = allelic_imbalence_vector[valid_samples]
        total_counts_vector = total_counts_vector[valid_samples]

        reg_site_indices = find_reg_site_positions_within_distance_cis(het_site_position, reg_site_positions, distance)

        # Now loop through all valid regulatory sites
        for reg_site_index in reg_site_indices:
            # Extract quantities of the regulatory site related to this iterations reg_site_index
            reg_site_position = reg_site_positions[reg_site_index]
            reg_site_id = reg_site_ids[reg_site_index]
            reg_site_vector = np.squeeze(np.asarray(reg_prob_mat[reg_site_index,:][:,valid_samples]))
            time_step_maf = get_heterozygous_maf(reg_site_vector)

            # Check to see if test is valid
            test_id = het_site_id + '_' + reg_site_id
            if test_id not in valid_tests:
                continue

            statistic, pvalue = run_statistical_test(reg_site_vector, allelic_imbalence_vector, statistical_test, normalization_method)
            #  Write test to ouput file
            t.write(str(chrom_num) + '\t' + het_site_id + '\t' + str(het_site_position) + '\t' + reg_site_id + '\t' + str(reg_site_position) + '\t' + str(time_step_maf) + '\t' + str(statistic) + '\t' + str(pvalue) + '\n')
    t.close()


##########################################
# Input Data
##########################################
het_prob_genotype_file = sys.argv[1]
allelic_counts_file = sys.argv[2]
output_file = sys.argv[3]
valid_tests_input_file = sys.argv[4]
min_reads = int(sys.argv[5])  # Minimum number of reads on heterozygous site to consider a sample valid
distance = float(sys.argv[6])  # Cis eqtl distance (EAGLE used 200 KB)
statistical_test = sys.argv[7]  # Statistical test to use in evaluating aseQTL (options are "wilcoxon" and "linear_regression")
normalization_method = sys.argv[8]  # Normalization method (only applicable if statistical_test="linear_regression", otherwise will be ignored). Options available are "standardize" and "quantile_normalize"
chrom_num = int(sys.argv[9])
time_step = int(sys.argv[10])

# We previously computed a list of variant-gene pairs that passed required filters
# Now extract dictionaries containing that information
valid_het_sites, valid_reg_sites, valid_tests = extract_valid_tests(valid_tests_input_file)

# Main driver function
driver(het_prob_genotype_file, allelic_counts_file, output_file, min_reads, distance, statistical_test, normalization_method, chrom_num, time_step, valid_het_sites, valid_reg_sites, valid_tests)
