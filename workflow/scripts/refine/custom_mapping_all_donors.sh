for smk in custom_mapping_*.smk; do snakemake -s $smk -j 50; done