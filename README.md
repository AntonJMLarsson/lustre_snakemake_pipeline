# Snakemake workflow: LUSTRE

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for processing LUSTRE data.

**Note (04/12/2024):** I'm currently refactoring the code to better comply with established best practices. Many components will change in the coming months.


## Usage

1. Run `demultiplex_fastq` to add the cell barcode to the fastq query name (see separarate repository). This requires the specification of the run folder, each lane can the be processed separately.

2. Concatenate the fastq files, and run `trim_and_map.smk`, which trims the reads using fastp and aligns the reads with bwa mem.

3. Run `preprocess.smk`, which does filtering, demultiplexing (separating reads per cell), sorts reads into groups based on KR, KNR and UNK.

4. Run `calling.smk`, which groups KNR and UNK reads into peaks and performs specific filtering steps to remove artefacts.

5. Run `reporting_regular.smk`, to produce a set of quality statistics before custom mapping.

6. Run `make_fastq_files.smk`, to identify partial reads which can be re-mapped for increased sensitivity.

7. For each donor in the dataset, run `make_custom_reference.smk` to obtain a unique reference to re-map to.

8. Run `custom_mapping_all_donors.smk`, to re-map for all donors separately.

9. Run `reporting_remapped.smk` to obtain the quality statistics after custom mapping.

10. Run `make_tables.smk` to make the insertion tables, which contain considerable information regarding each detected insertion, and read count tables.

# TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.
* Describe other auxillary pipelines contained in this repository.
* Describe insertion table output.
