#!/bin/bash
#SBATCH --time=05:00:00 --mem=10GB

dosage_genotype_file="$1"
quantile_normalized_expression="$2"
gencode_gene_annotation_file="$3"
output_dir="$4"
eqtl_distance="$5"
maf_cutoff="$6"
time_step="$7"
visualize_independent_time_step_eqtl_dir="$8"
full_ipsc_qtl_data="$9"
full_heart_eqtl_data="${10}"

output_file_prefix=$output_dir"eqtl_prepare_eqtl_distance_"$eqtl_distance"_maf_cutoff_"$maf_cutoff

if false; then
# Perorm qtl analysis in each chromosome seperately (okay because doing cis eqtls)
for chrom_num in $(seq 1 22); do
    output_root=$output_file_prefix"_time_step_"$time_step"_chrom_"$chrom_num"_"
    
    # Prepare input files to matrix eqtl
    python prepare_independent_time_step_matrix_eqtl_files.py $time_step $chrom_num $dosage_genotype_file $quantile_normalized_expression $gencode_gene_annotation_file $output_root $maf_cutoff
    # Input file names
    matrix_eqtl_genotype_file=$output_root"genotype.txt"
    matrix_eqtl_variant_loc_file=$output_root"variant_location.txt"
    matrix_eqtl_expression_file=$output_root"expression.txt"
    matrix_eqtl_gene_location_file=$output_root"gene_location.txt"
    matrix_eqtl_output_file=$output_root"matrix_eqtl_out.txt"

    # Run matrix eqtl
    Rscript matrix_eqtl_wrapper.R $matrix_eqtl_genotype_file $matrix_eqtl_variant_loc_file $matrix_eqtl_expression_file $matrix_eqtl_gene_location_file $eqtl_distance $matrix_eqtl_output_file
    # Convert matrix eqtl output to easier to use format
    python convert_matrix_eqtl_output_to_standard_association_output.py $matrix_eqtl_output_file $output_root"eqtl_results.txt" $matrix_eqtl_genotype_file $matrix_eqtl_variant_loc_file $matrix_eqtl_gene_location_file $chrom_num
    
    # Delete temporary files
    rm $matrix_eqtl_genotype_file
    rm $matrix_eqtl_variant_loc_file
    rm $matrix_eqtl_expression_file
    rm $matrix_eqtl_gene_location_file
    rm $matrix_eqtl_output_file
done
# Multiple testing correction for genome wide significance
python independent_time_step_multiple_testing_correction.py $output_file_prefix"_time_step_"$time_step"_" "eqtl_results.txt"



# Visualize distribution of pvalues from eqtl analyses at each of the time steps independently
Rscript visualize_eqtls.R $output_file_prefix"_time_step_"$time_step"_bonferonni_correction.txt" $output_file_prefix"_time_step_"$time_step"_eqtl_results.txt" $visualize_independent_time_step_eqtl_dir"eqtl_distance_"$eqtl_distance"_maf_cutoff_"$maf_cutoff"_time_step_" $time_step
fi

# We are given significant eqtls from Nick Banovich's ipsc eqtl analysis on about 60 subjects.
# We extract both:
#  1. the p-values of those eqtls in our data.
#  2. p-values of matched snp-gene pairs (matched for distance to tss and maf)
version="ipsc"
python get_snp_gene_pairs_from_other_data_and_match_for_background.py $output_file_prefix"_time_step_"$time_step"_eqtl_results.txt" $visualize_independent_time_step_eqtl_dir"ipsc_banovich_comparison_distance_"$eqtl_distance"_maf_cutoff_"$maf_cutoff"_time_step_"$time_step"_" $full_ipsc_qtl_data $eqtl_distance $version

if false; then
# Do same thing as above, except for GTEx v7 heart left ventricle eqtls
version="heart"
python get_snp_gene_pairs_from_other_data_and_match_for_background.py $output_file_prefix"_time_step_"$time_step"_eqtl_results.txt" $visualize_independent_time_step_eqtl_dir"gtex_v7_heart_left_ventricle_comparison_distance_"$eqtl_distance"_maf_cutoff_"$maf_cutoff"_time_step_"$time_step"_" $full_heart_eqtl_data $eqtl_distance $version
fi