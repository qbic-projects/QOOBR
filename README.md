# Data analysis scripts for QMKMK - CLAD project

This repository contains a conda environment file to reproduce all the data analysis results under `conda_environment.yml`.

## BCR analysis

BCR analysis was performed with the nf-core/airrflow pipeline v2.0.0 (https://github.com/nf-core/airrflow).

```bash
nextflow run nf-core/bcellmagic -r 2.0.0 -profile cfc --input samplesheet.tsv --protocol pcr_umi --cprimer_position R1 --umi_length 8 --vprimers VPrimers.fasta --cprimers CPrimers_IG.fasta --index_file --singularity_pull_docker_containers --set_cluster_threshold --cluster_threshold 0.079
```

The repertoire analysis results from the pipeline can be found under:

- BCR_analysis/repertoire_per_patient: Pipeline run defining the clones per patient. Using as samplesheet: `repertoire_per_patient/samplesheet.tsv`
- BCR_analysis/repertoire_per_time_point: Pipeline run defining the clones per time point. Using as samplesheet: `repertoire_per_timepoint/samplesheet.tsv`

The report directly outputted by the pipeline was enhanced to have additional plots. The enhanced report can be found under `repertoire_per_timepoint/repertoire_comparison_timepoint.Rmd` and `repertoire_per_timepoint/repertoire_comparison_timepoint.nb.html`, containing all figures used for the publication.

### Clonal overlap analysis

The scripts for clonal overlap analysis for exclusive and inclusive intersects can be found under `repertoire_per_patient/Clone_overlaps_bootstrapped_exclusive_intersects.R` and `repertoire_per_patient/Clone_overlaps_bootstrapped.R`, respectively.

## Proteomics analysis

Proteomics data was analyzed with MaxQuant, providing as reference database the nucleotide sequences of each of the subjects. The output of MaxQuant was further processed with custom scripts.

- `Proteomics_analysis/translation/` contains the scripts used for translating the BCR nucleotide sequences into amino acid sequences.
- `peptide_matches_CDR3.R`: filtered the peptide matches to remove multimatches and filter matches to CDR regions.
- `proteome-transcriptome-overlaps.Rmd`: filtering peptide matches according to CDR regions matches and plots for figures.
