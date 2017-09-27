import numpy as np
import os
import sys
import pdb

# Method to run aseqtl analysis on one chromosome
#INPUTS:
## 1. output_file_prefix: Prefix to outputfile
## 2. allelic_imbalence_mat: matrix of dim num_het_sites X samples. Has na in elements where het_site, sample are not heterozygous or have 0 reads mapping
## 3. total_counts_mat: matrix of dim num_het_sites X samples. Has na in elements where het_site, sample are not heterozygous
## 4. het_site_positions: list of integers of length num_het_sites where each element is the position of the corresponding heterozygous site
## 5. het_site_ids: List of strings of length num_het_sites where each element is the name of the corresponding heterozygous site
## 6. reg_prob_matrix: matrix of dim num_regulatory_sites X samples. Where each element is the heterozygous probability that the site, sample is heterozygous
## 7. reg_site_positions: List of integers of length num_regulatory_sites where each element is the positions of the corresponding regulatory site
## 8. reg_site_ids: List of strings of length num_regulatory_sites where each element is the name of the corresponding regulatory site
## 9. min_reads: Minimum number of total reads a (het_site,sample) must have to be considered for this analysis
## 10. min_samples_per_time_step: Minimum number of samples that must pass filters for a specific test
## 11. min_fraction_in_test_group: Minimum fraction of valid samples with less popular version of regulatory variant (homozygous reg variant vs heterozygous reg variant)
## 12. distance: cis-regulatory distance
## 13. statistical_test: Statistical test to use in evaluating aseQTL (options are "wilcoxon" and "linear_regression")
## 14. normalization_method: Normalization method (only applicable if statistical_test="linear_regression", otherwise will be ignored). Options available are "standardize" and "quantile_normalize"
def run_ase_qtl(output_file_prefix, allelic_imbalence_mat, total_counts_mat, het_site_positions, het_site_ids, reg_prob_matrix, reg_site_positions, reg_site_ids, min_reads, min_samples_per_time_step, min_fraction_in_test_group, distance, statistical_test, normalization_method):
    pdb.set_trace()