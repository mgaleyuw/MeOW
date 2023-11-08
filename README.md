# Methylation Operation Wizard (MeOW) 0.2.3

A fast identifier of differentially methylated regions of 5mC called BAM files. The magnitude and statistical significance of difference is assessed using either corrected t-tests or beta regression and results are visualized in R Notebooks.

# Basic Setup

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

In this version, the wildcard SAMPLEID is expected as part of the bam file name and is constrained to match the pattern `M\d\d\d\d`. This pattern can be changed in `Snakefile`, but the wildcard must be constrained.

**MeOW will run slowly the first time as conda packages are downloaded and installed.**

# Using Pre-existing Control Datasets

To use exisiting control databases, copy the contents of `optional_pregenerated_results/{5mcmodel}/results` to `results` and `optional_pregenerated_results/{5mcmodel}/config` to `config/` and `workflow/config`

# Results

The following files are produced when running t-tests or beta regressions:
1. `results/meow.reference.annotated.tsv`: This file is generated the first time MeOW is run using the controls listed in your metadata. It contains the mean methylation frequencies at each CpG position within the intervals defined in your config file. The default is to calculate these within the Hg38 CpG islands defined by UCSC.
2. `results/meow.reference.index.tsv` : This is an index file
3. `results/{SAMPLEID}.cpg.methyl.indexed.tsv` : This is a two column file with methylation frequencies for the sample SAMPLEID and indexed positions.
- `results/{SAMPLEID}_meow_results`
    1. `{SAMPLEID}.{testype}.dmr.html` : This is an HTML report that summarizes statistically significant test results. If there are no statistically significant results, this file may fail to generate. [This should be fixed]
    2.  `{SAMPLEID}.{testtype}.dmr.Rmd` : This is an Rnotebook that is always generated regardless of whether or not test results are significant. It can be run in RStudio and edited to widen parameters in the case of no significant results.
    3.  `{SAMPLEID}_all_DMR_{testtype}.csv` : All bootstrapped test results for regions (CpG islands) included in analysis. Includes effect size (Cohen's d), corrected test statistic, gene annotation, and CpG label.
    4. `{SAMPLEID}_all_DMR_{testtype}.csv` : Results of DMR test filtered to p less than or equal to 0.05 and effect magnitudes of 1.75.
    5. `{SAMPLEID}_all_DMR_{testtype}_ribbons.png` : Plot of differentially methylated regions in test sample compared to the confidence interval around control samples for the same positions.
    6. `{SAMPLEID}_all_DMR_{testtype}_boxplots.png` : As above but with boxplots.
    7. `{SAMPLEID}_all_DMR_{testtype}_biplot.png` : A biplot of all results, where comparisons that are both statistically significant and with large effect sizees are highlighted in red and labeled.

# Conda

MeOW can be used without conda (don't include `--use-conda` when running Snakemake), as long as the following requirements are met:
- python >= 3.7
- snakemake
- samtools=1.17
- rust=1.72.1
- r-base=4.2.1
    - r-essentials
    - r-ggplot2
    - r-plyr
    - r-dplyr
    - r-scales
    - r-stringr
    - r-cowplot
    - r-tidyr
    - r-reshape2
    - r-fitdistrplus
    - r-base64enc
    - r-digest
    - r-evaluate
    - r-glue
    - r-highr
    - r-htmltools
    - r-jsonlite
    - r-knitr
    - r-magrittr
    - r-markdown
    - r-mime
    - r-stringi
    - r-tinytex
    - r-xfun
    - r-yaml
    - r-ggrepel