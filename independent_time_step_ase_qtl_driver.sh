#!/bin/bash
#SBATCH --time=10:00:00 --mem=10G

het_prob_genotype_file="$1"
allelic_counts_file="$2"
output_stem="$3"
min_reads="$4"
min_samples_per_time_step="$5"
min_fraction_in_test_group="$6"
distance="$7"
statistical_test="$8"
normalization_method="$9"
min_fraction_biallelic="${10}"
visualize_independent_time_step_ase_qtl_dir="${11}"
qtl_dir="${12}"



output_qtl_file_prefix=$qtl_dir$output_stem
output_visualization_file_prefix=$visualize_independent_time_step_qtl_dir$output_stem

#######################################################################################
# The goal of this script is to create a list of variant-gene pairs that pass all
# of our ase-qtl filters.
# The filters are as follows:
# A variant gene pair must pass all of the following in EVERY time step:
#### 1. distance between reg. variant and heterozygous site is < distance
#### 2. There are at least min_samples_per_time_step that are heterozygous at the exonic variant
#### 3. Of the the heterozygous samples, at least min_fraction_biallelic % of them show biallelic expression
####- where biallelic expression is both: a) greater than min_reads total reads and b) allelic imbalence < .49
#### 4. Regulatory variant must have MAF >= min_fraction_in_test_group
########################################################################################
# Do this for each chromosome seperately
if false; then
for chrom_num in $(seq 1 22); do
    python extract_valid_ase_tests_across_time_steps.py $het_prob_genotype_file $allelic_counts_file $output_qtl_file_prefix"chrom_"$chrom_num"_valid_tests.txt" $min_reads $min_samples_per_time_step $min_fraction_in_test_group $distance $chrom_num $min_fraction_biallelic
done
fi


for time_step in $(seq 0 15); do
    if false; then
    for chrom_num in $(seq 1 22); do
        echo "Time step: "$time_step  "Chrom_num: "$chrom_num
        python independent_time_step_ase_qtl_chromosome_specific.py $het_prob_genotype_file $allelic_counts_file $output_qtl_file_prefix"time_step_"$time_step"_chrom_"$chrom_num"_aseqtl_results.txt" $output_qtl_file_prefix"chrom_"$chrom_num"_valid_tests.txt" $min_reads $distance $statistical_test $normalization_method $chrom_num $time_step
    done
    python independent_time_step_multiple_testing_correction.py $output_qtl_file_prefix"time_step_"$time_step"_" "aseqtl_results.txt"
    fi
    # Visualize distribution of pvalues from eqtl analyses at each of the time steps independently
    Rscript visualize_eqtls.R $output_qtl_file_prefix"time_step_"$time_step"_bonferonni_correction.txt" $output_qtl_file_prefix"time_step_"$time_step"_aseqtl_results.txt" $visualize_independent_time_step_ase_qtl_dir$output_stem $time_step
done
