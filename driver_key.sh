
# This scripts assumes you have run the ipsc_preproccess_pipeline first (https://github.com/BennyStrobes/ipsc_preprocess_pipeline)
# And have save the results here:
preprocess_dir="/project2/gilad/bstrober/ipsc_differentiation/preprocess/"


##########################################################################
# Input Files
##########################################################################

# Heterozygous probability genotype file for all cell_lines in our analysis
# These heterozygous probabilities come from impute2
het_prob_genotype_file=$preprocess_dir"genotype/YRI_het_prob_genotype.vcf"

# Dosage genotype for all cell lines in our analysis
dosage_genotype_file=$preprocess_dir"genotype/YRI_genotype.vcf"

# Processed Allelic counts
# Assuming heterozygous threshold of .999 in terms of calling heterozygous sites
allelic_counts_file=$preprocess_dir"processed_allelic_counts/allelic_counts_gene_mapped_het_prob_999.txt"

# SVA loadings for each of the samples
# This is our covariate file
sva_loading_file=$preprocess_dir"covariates/sva_loadings.txt"


##########################################################################
# Output Directories (Assumes these directories exist before script starts)
##########################################################################
# Root directory for ipsc qtl pipelines
ipsc_qtl_root="/project2/gilad/bstrober/ipsc_differentiation/qtl_pipelines/"

# Output directory for running ase qtl analysis at each time step independently
independent_time_step_qtl_dir=$ipsc_qtl_root"independent_time_step_qtl/"






##########################################################################
# Output Files (Do not exist before program starts). 
##########################################################################





##########################################################################
# Scripts/Functions/Analyses
##########################################################################

############ Run independent time step qtl analysis
##MODEL PARAMETERS
# Cis eqtl distance (EAGLE used 200 KB)
distance="200000"
# Minimum number of reads on heterozygous site to consider a sample valid
min_reads="5"
# Minimum number of valid samples we will run a qtl test on
min_samples_per_time_step="7"
# Minimum fraction of valid samples with less popular version of regulatory variant (homozygous reg variant vs heterozygous reg variant)
min_fraction_in_test_group=".2"
# Statistical test to use in evaluating aseQTL
# Options available are "wilcoxon" and "linear_regression
statistical_test="wilcoxon"
# Normalization method (only applicable if statistical_test="linear_regression", otherwise will be ignored)
# Options available are "standardize" and "quantile_normalize"
normalization_method="na"

## Output file prefix
output_file_prefix=$independent_time_step_qtl_dir"independent_time_step_qtl_distance_"$distance"_min_reads_"$min_reads"_min_samples_"$min_samples_per_time_step"_min_fraction_in_test_group_"$min_fraction_in_test_group"_"$statistical_test"_normalization_method_"$normalization_method"_"

sh independent_time_step_qtl_driver.sh $het_prob_genotype_file $allelic_counts_file $output_file_prefix $min_reads $min_samples_per_time_step $min_fraction_in_test_group $distance $statistical_test $normalization_method





