import numpy as np
import os
import sys
import pdb
from scipy.stats import ranksums
from scipy.stats import rankdata
from scipy import stats
import statsmodels.api as sm




# Find indices of reg_site_positions that are within distance of het_site_position
def find_reg_site_positions_within_distance(het_site_position, reg_site_positions, distance):
    # Compute absolute distance from each reg_site and the current heterozygous site
    absolute_difference_vector = abs(reg_site_positions - het_site_position)
    # Find all indices where distance between regulatory site and heterozygous site is less than distance
    valid_indices = np.where(absolute_difference_vector <= distance)[0]
    return valid_indices

# Make sure >= min_fraction_in_test_group of the samples belong to either homozygous reg variant or heterozygous reg variant
def pass_min_fraction_in_test_group_filter(reg_site_vector, min_fraction_in_test_group, statistical_test):
    # Wilcoxon hard codes assignments
    if statistical_test == 'wilcoxon':
        num_heterozygous = float(len(np.where(reg_site_vector >= .5)[0]))
        num_homozygous = float(len(reg_site_vector) - num_heterozygous)
    # linear regression uses heterozygous probabilities
    if statistical_test == 'linear_regression':
        num_heterozygous = np.sum(reg_site_vector)
        num_homozygous = len(reg_site_vector) - num_heterozygous
    # Initialize boolean filter pass to true
    pass_filter = True
    total_samples = num_heterozygous + num_homozygous
    # check if it does not pass filter
    if num_heterozygous/total_samples < min_fraction_in_test_group or num_homozygous/total_samples < min_fraction_in_test_group:
        pass_filter = False
    return pass_filter

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


# Method to run aseqtl analysis on one chromosome
#INPUTS:
## 1. output_file_prefix: Prefix to outputfile
## 2. allelic_imbalence_mat: matrix of dim num_het_sites X samples. Has na in elements where het_site, sample are not heterozygous or have 0 reads mapping
## 3. total_counts_mat: matrix of dim num_het_sites X samples. Has na in elements where het_site, sample are not heterozygous
## 4. het_site_positions: list of integers of length num_het_sites where each element is the position of the corresponding heterozygous site
## 5. het_site_ids: List of strings of length num_het_sites where each element is the name of the corresponding heterozygous site
## 6. reg_prob_mat: matrix of dim num_regulatory_sites X samples. Where each element is the heterozygous probability that the site, sample is heterozygous
## 7. reg_site_positions: List of integers of length num_regulatory_sites where each element is the positions of the corresponding regulatory site
## 8. reg_site_ids: List of strings of length num_regulatory_sites where each element is the name of the corresponding regulatory site
## 9. min_reads: Minimum number of total reads a (het_site,sample) must have to be considered for this analysis
## 10. min_samples_per_time_step: Minimum number of samples that must pass filters for a specific test
## 11. min_fraction_in_test_group: Minimum fraction of valid samples with less popular version of regulatory variant (homozygous reg variant vs heterozygous reg variant)
## 12. distance: cis-regulatory distance
## 13. statistical_test: Statistical test to use in evaluating aseQTL (options are "wilcoxon" and "linear_regression")
## 14. normalization_method: Normalization method (only applicable if statistical_test="linear_regression", otherwise will be ignored). Options available are "standardize" and "quantile_normalize"
## 15. chrom_num: integer corresponding to which chromosome the analysis was done on 
## 16. min_fraction_biallelic: Minimum fraction of heterozygous samples that show biallelic expression
def run_ase_qtl(output_file_prefix, allelic_imbalence_mat, total_counts_mat, het_site_positions, het_site_ids, reg_prob_mat, reg_site_positions, reg_site_ids, min_reads, min_samples_per_time_step, min_fraction_in_test_group, distance, statistical_test, normalization_method, chrom_num, min_fraction_biallelic):
    # Get arrays and matrices into correct format
    reg_site_positions = np.asarray(reg_site_positions)
    allelic_imbalence_mat = np.asmatrix(allelic_imbalence_mat)
    total_counts_mat = np.asmatrix(total_counts_mat)
    reg_prob_mat = np.asmatrix(reg_prob_mat)

    # Some simple checks
    num_het_sites, num_samples = allelic_imbalence_mat.shape
    num_reg_sites, num_samples2 = reg_prob_mat.shape
    if num_samples != num_samples2:
        print('FATAL ERROR: Matrix shape mismatch between alllelic_imbalence_mat and reg_prob_matrix')
    
    # Open output file handle
    t = open(output_file_prefix + '_aseqtl_results.txt', 'w')
    # Write heaer
    t.write('chom_num\thet_site_id\thet_site_position\tregulatory_site_id\tregulatory_site_position\ttest_statistic\tpvalue\n')

    # First loop through all instances of the dependent variable (allelic imbalence)
    for het_site_index in range(num_het_sites):
        # Extract quantities of allelic_imbalence related to this iterations het_site_index
        allelic_imbalence_vector = np.squeeze(np.asarray(allelic_imbalence_mat[het_site_index,:]))
        total_counts_vector = np.squeeze(np.asarray(total_counts_mat[het_site_index,:]))
        het_site_position = het_site_positions[het_site_index]
        het_site_id = het_site_ids[het_site_index]

        ##FILTERS
        # Get number of samples that have a heterozygous site
        num_samples_with_het_site = np.sum(np.isnan(total_counts_vector)==False)
        # Make sure there are at least min_samples_per_time_step valid_samples
        if num_samples_with_het_site < min_samples_per_time_step:
            continue
        # Extract indices of samples that have at least min_reads reads overlapping site and do not do not show allelic imbalence
        valid_samples = (total_counts_vector >= min_reads)*(allelic_imbalence_vector < .49)
        num_valid_samples = np.sum(valid_samples)
        # Make sure the fraction of biallelic samples (out of the number of heterzygous samples) is >= min_fraction_biallelic
        if num_valid_samples/float(num_samples_with_het_site) < min_fraction_biallelic:
            continue

        allelic_imbalence_vector_filtered_samples = allelic_imbalence_vector[valid_samples]
        # Find indices of reg_site_positions that are within distance of het_site_position
        reg_site_indices = find_reg_site_positions_within_distance(het_site_position, reg_site_positions, distance)
        # Now loop through all valid regulatory sites
        for reg_site_index in reg_site_indices:
            # Extract quantities of the regulatory site related to this iterations reg_site_index
            reg_site_position = reg_site_positions[reg_site_index]
            reg_site_id = reg_site_ids[reg_site_index]
            reg_site_vector = np.squeeze(np.asarray(reg_prob_mat[reg_site_index,:][:,valid_samples]))
            # Make sure >= min_fraction_in_test_group of the samples belong to either homozygous reg variant or heterozygous reg variant
            if pass_min_fraction_in_test_group_filter(reg_site_vector, min_fraction_in_test_group, statistical_test) == False:
                continue
            # Add this point the allelic imbalence (at het. snp) and the regulatory variant have passed all of our filters. We can now run analysis
            statistic, pvalue = run_statistical_test(reg_site_vector, allelic_imbalence_vector_filtered_samples, statistical_test, normalization_method)
            #  Write test to ouput file
            t.write(str(chrom_num) + '\t' + het_site_id + '\t' + str(het_site_position) + '\t' + reg_site_id + '\t' + str(reg_site_position) + '\t' + str(statistic) + '\t' + str(pvalue) + '\n')
    # Close output handle
    t.close()