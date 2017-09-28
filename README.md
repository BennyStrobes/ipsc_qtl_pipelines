# ipsc qtl pipelines

This pipeline runs multiple ase-qtl analysis on the ipsc allelic count data. All analyses can be controlled through changing 'if false; then' commands in 'driver_key.sh' The types of ase-qtl analysis run are:

	1. independent time step aseq-qtl: Run standard ase qtl analysis for each time step independently. Implemented using wilcoxon or linear regression.


## Dependencies

These scripts assumes you have already preprocessed the rna-seq data using [ipsc_preprocess_pipeline](https://github.com/BennyStrobes/ipsc_preprocess_pipeline).

## Deliverables

As for important output files from this pipeline (using dir_names defined in 'driver_key.sh':

	1. $independent_time_step_qtl_dir"independent_time_step_qtl_distance_"$distance"_min_reads_"$min_reads"_min_fraction_biallelic_"$min_fraction_biallelic"_min_samples_"$min_samples_per_time_step"_min_fraction_in_test_group_"$min_fraction_in_test_group"_"$statistical_test"_normalization_method_"$normalization_method"_"

The variables in the output file refer to parameters in the model. See `driver_key.sh` for a description of all of the parameters.


Let me know if something isn't clear.

## Computer cluster

This pipeline was written to run on midway rcc

## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes)