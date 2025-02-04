#!/usr/bin/bash

# replace COHORT with desired cohort to run

COHORT=Mexico

snakemake \
  --use-conda \
  --profile ./resources/profiles/smk-simple-slurm/simple 
  --snakefile Snakefile-bin \
  --configfile config.${COHORT}.yaml \
  --jobs 100
