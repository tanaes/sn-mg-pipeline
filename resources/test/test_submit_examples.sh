# test local submission
snakemake -j 8 \
    --local-cores 8 \
    --use-conda \
    --conda-prefix /home/jgsanders/.sn-envs

# test local submission, no taxonomy profiling
snakemake -j 8 \
    --local-cores 8 \
    --use-conda \
    --conda-prefix /home/jgsanders/.sn-envs \
    no_profile

# test SLURM submission, no taxonomy profiling
snakemake \
    --profile resources/profiles/smk-simple-slurm/simple \
    --local-cores 8 \
    --use-conda \
    --conda-prefix /home/jgsanders/.sn-envs \
    --cluster-status resources/profiles/smk-simple-slurm/extras/status-sacct.sh \
    no_profile