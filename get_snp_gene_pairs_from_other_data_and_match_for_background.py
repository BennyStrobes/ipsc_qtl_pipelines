import numpy as np
import os
import sys
import pdb
import random
import math

# Return the bin number corresponding to this distance
def get_distance_bin(distance, distance_bin_size):
    return int(math.floor(distance/distance_bin_size))


# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))

# Make object that input you can input distance and maf into, and it will return an array of all the pvalues contained in that bin
def make_background_qtls_object(our_eqtl_file, distance_bin_size, maf_bin_size, eqtl_distance):
    ####################
    # Initialize object
    ####################
    background_qtls = []
    # number of bins needed for maf and distance
    num_distance_bins = int(math.ceil(eqtl_distance/distance_bin_size + 1))
    num_maf_bins = int(math.ceil(.5/maf_bin_size + 1))
    # Add each possible bin
    for distance_bin in range(num_distance_bins):
        background_qtls.append([])
        for maf_bin in range(num_maf_bins):
            background_qtls[distance_bin].append([])

    ####################
    # Add pvalues from our_eqtl_file to the appropriate bin
    ####################
    f = open(our_eqtl_file)
    head_count = 0  # Used to skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1  # Skip Header
            continue
        # Extract relevent fields from column delimited file
        gene_position = float(data[2])
        variant_position = float(data[4])
        distance = abs(gene_position-variant_position)
        maf = float(data[5])
        pvalue = float(data[7])
        rs_id = data[3]
        if maf > .5 or maf < .1:
            print('maf errro')
            pdb.set_trace()

        # Check to make sure we are only dealing with labeled rs_ids
        # This should have been taken care of in `prepare_independent_time_step_matrix_eqtl_files.py`. But just checking here
        if rs_id == '.':
            print('ERROR')
            continue

        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)

        # Add pvalue of this variant gene pair to the appropriate bin
        background_qtls[distance_bin][maf_bin].append(pvalue)
    return background_qtls

# Extract dictionary list of variant_geneID pairs that are found to be eqtls in reference_eqtl_file
def extract_reference_eqtls(reference_eqtl_file,version):
    eqtls = {}
    # Stream file
    f = open(reference_eqtl_file)
    head_count = 0  # used to skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:  # Skip header
            head_count = head_count + 1
            continue
        # Extract gene id and variant id in different ways depending on where the data comes from
        if version == 'ipsc':
            # Extract relevent fields from line
            gene_id = data[0].split('.')[0]
            rsid = data[1].split('.')[0]
        elif version == 'heart':
            gene_id = data[0].split('.')[0]
            rsid = data[18]
        # Concatenate geneID and rsID to get name of test
        test_name = gene_id + '_' + rsid
        # Add eqtl to dictionary
        eqtls[test_name] = 1
    return eqtls

# Stream our_eqtl_file and find those that are in our_eqtl_file
# For each one of those hits, randomly select background pvalue
def extract_and_print_pvalues(background_qtls, reference_eqtls, our_eqtl_file, output_file):
    # Initialize output handle
    t = open(output_file, 'w')
    t.write('real_pvalue\tmatched_pvalue\n')

    # Stream our_eqtl_file
    head_count = 0
    f = open(our_eqtl_file)
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1  # Skip Header
            continue
        # Extract relevent fields
        pvalue = float(data[7])
        gene_id = data[1]
        rs_id = data[3]
        test_id = gene_id + '_' + rs_id
        gene_position = float(data[2])
        variant_position = float(data[4])
        distance = abs(gene_position-variant_position)
        maf = float(data[5])
        # Check to see if test_id in reference_eqtls
        if test_id not in reference_eqtls:
            continue
        # Now randomly select a background variant-gene pair's pvalue for reference
        # Return the bin number corresponding to this distance
        distance_bin = get_distance_bin(distance, distance_bin_size)
        # Return the bin number corresponding to this distance
        maf_bin = get_maf_bin(maf, maf_bin_size)
        
        # Extract list of pvalue corresponding to this distance bin and maf_bin
        pvalues = background_qtls[distance_bin][maf_bin]

        if len(pvalues) < 5:
            print('small bin error')
            print(len(pvalues))
            print(distance_bin)
            print(maf_bin)
            print(distance)
            print(maf)
        # Randomly select one element from the list
        random_pvalue = random.choice(pvalues)

        # Write to output file
        t.write(str(pvalue) + '\t' + str(random_pvalue) + '\n')




###############################################################################################
# Objective
###############################################################################################
# We are given significant eqtls from another eqtl analysis.
# We extract both:
#  1. the p-values of those eqtls in our data.
#  2. p-values of matched snp-gene pairs (matched for distance to tss and maf)



###############################################################################################
# Input Data
###############################################################################################
# eqtl results from our analysis
our_eqtl_file = sys.argv[1]
# Prefix to output files
output_file_root = sys.argv[2]
# eqtl results we are treating as gold standard
reference_eqtl_file = sys.argv[3]
# Max distance a variant and gene can be from one another
eqtl_distance = float(sys.argv[4])
# Whether this analysis is being run for:
## 1. 'ipsc' : nick banovich's file format
## 2. 'heart' : gtex v7 file format
version = sys.argv[5]


# Size of bins used for background matching
distance_bin_size = 14000
maf_bin_size = .08


# Make object that input you can input distance and maf into, and it will return an array of all the pvalues contained in that bin
background_qtls = make_background_qtls_object(our_eqtl_file, distance_bin_size, maf_bin_size, eqtl_distance)


# Extract dictionary list of geneID_variantID pairs that are found to be eqtls in reference_eqtl_file
reference_eqtls = extract_reference_eqtls(reference_eqtl_file,version)

# Output_file to save results to
output_file = output_file_root + 'real_v_matched_controls.txt'

# Stream our_eqtl_file and find those that are in our_eqtl_file
# For each one of those hits, randomly select background pvalue
extract_and_print_pvalues(background_qtls, reference_eqtls, our_eqtl_file, output_file)