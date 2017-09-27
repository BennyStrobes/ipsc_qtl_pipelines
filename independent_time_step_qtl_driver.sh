
het_prob_genotype_file="$1"
allelic_counts_file="$2"
output_file_prefix="$3"
min_reads="$4"
min_samples_per_time_step="$5"
min_fraction_in_test_group="$6"
distance="$7"
statistical_test="$8"
normalization_method="$9"


chrom_num="1"
time_step="2"

python independent_time_step_qtl_chromosome_specific.py $het_prob_genotype_file $allelic_counts_file $output_file_prefix$time_step"_"$chrom_num"_" $min_reads $min_samples_per_time_step $min_fraction_in_test_group $distance $statistical_test $normalization_method $chrom_num $time_step