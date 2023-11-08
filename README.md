# Snakemake Version of MeOW

Miranda Perez Galey Zalusky

To use, set up configuration like so:

`/config` should contain:
 -  A tab separated file with sample IDs, the full path to an aligned BAM file with MM/ML tags encoding methylation, and a designation as Test/Control. An example file is included (`config/blood_5mc_bam_files_20231023.tsv`)
 - A comma separated file with minimum metadata for control files (example included, `config/blood_5mc_control_metadata_20231023.tsv`)
 - A comma separated file with minimum metadata for test samples (example included, `config/test_metadata.csv`)
 - `targets.txt`, listing the sample IDs for which you wish to generate meow tests, assuming you wish to use the wild card based target rules.
 - `config.yaml` : update this file with appropriate filenames of the above and choose the methylation profile you would like to use. Controls and Test samples must use the same methylation profiles.

 Place additional copies of the metadata files in `workflow/config/`

Within a minimal smakemake conda environment: run from this directory using `snakemake --use-conda --cores N` where N=the number of threads you would like to use. It is reccommended to run MeOW with at least 25 cores.

Results will be generated for each test sample in /results/{SAMPLEID}_meow_results/

