#!/bin/bash
#SBATCH --time=10:00:00 --mem=20G

het_prob_genotype_file="$1"
allelic_counts_file="$2"
output_file_prefix="$3"
min_reads="$4"
min_samples_per_time_step="$5"
min_fraction_in_test_group="$6"
distance="$7"
statistical_test="$8"
normalization_method="$9"
min_fraction_biallelic="${10}"




for time_step in $(seq 0 15); do
    for chrom_num in $(seq 1 22); do
        echo "Time step: "$time_step  "Chrom_num: "$chrom_num
        python independent_time_step_qtl_chromosome_specific.py $het_prob_genotype_file $allelic_counts_file $output_file_prefix"time_step_"$time_step"_chrom_"$chrom_num $min_reads $min_samples_per_time_step $min_fraction_in_test_group $distance $statistical_test $normalization_method $chrom_num $time_step $min_fraction_biallelic
    done
    python independent_time_step_multiple_testing_correction.py $output_file_prefix"time_step_"$time_step"_"
done