from os.path import basename, dirname, join
from shutil import copyfile
from glob import glob

localrules: consolidate_DAS_Tool_bins, prepare_dRep

rule metabat2_Fasta_to_Contig2Bin:
    """
    Uses Fasta_to_Contig2Bin script in DAS Tools to generate a contigs2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/metabat2/{mapper}/run_metabat2/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        contigs2bin="output/selected_bins/metabat2/{mapper}/contigs2bin/{contig_sample}_contigs2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/metabat2/{mapper}/contigs2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/metabat2/{mapper}/contigs2bin/{contig_sample}.log"
    shell:
        """ 
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fa > {output.contigs2bin} 2> {log}
        """


rule maxbin2_Fasta_to_Contig2Bin:
    """
    Uses Fasta_to_Contig2Bin script in DAS Tools to generate a contigs2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/maxbin2/{mapper}/run_maxbin2/{contig_sample}/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        contigs2bin="output/selected_bins/maxbin2/{mapper}/contigs2bin/{contig_sample}_contigs2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/maxbin2/{mapper}/contigs2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/maxbin2/{mapper}/contigs2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fasta > {output.contigs2bin} 2> {log}
        """

rule concoct_Fasta_to_Contig2Bin:
    """
    Uses Fasta_to_Contig2Bin script in DAS Tools to generate a contigs2bin.tsv file.
    """
    input:
        bins = lambda wildcards: expand("output/binning/concoct/{mapper}/extract_fasta_bins/{contig_sample}_bins/",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample)
    output:
        contigs2bin="output/selected_bins/concoct/{mapper}/contigs2bin/{contig_sample}_contigs2bin.tsv"
    conda:
        "../env/selected_bins.yaml"
    benchmark:
        "output/benchmarks/selected_bins/concoct/{mapper}/contigs2bin/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/concoct/{mapper}/contigs2bin/{contig_sample}.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fa  | perl -pe 's/(.*?)\s.*\\t(.*)$/\\1\\t\\2/' > {output.contigs2bin} 2> {log}
        """

rule run_DAS_Tool:
    """
    Selects bins using DAS_Tool
    """
    input:
        metabat2 = lambda wildcards: expand("output/selected_bins/metabat2/{mapper}/contigs2bin/{contig_sample}_contigs2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        maxbin2 = lambda wildcards: expand("output/selected_bins/maxbin2/{mapper}/contigs2bin/{contig_sample}_contigs2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        concoct = lambda wildcards: expand("output/selected_bins/concoct/{mapper}/contigs2bin/{contig_sample}_contigs2bin.tsv",
                mapper = config['mappers'],
                contig_sample = wildcards.contig_sample),
        contigs = lambda wildcards: expand("output/assemble/{assembler}/{contig_sample}.contigs.fasta",
                    assembler = config['assemblers'],
                    contig_sample = wildcards.contig_sample)
    output:
        out="output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}/{contig_sample}_DASTool_summary.tsv"
    params:
        basename = "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}/{contig_sample}",
        search_engine = config['params']['das_tool']['search_engine']
    conda:
        "../env/selected_bins.yaml"
    threads:
        res['run_DAS_Tool']['threads']
    resources:
        partition = res['run_DAS_Tool']['partition'],
        mem_mb = res['run_DAS_Tool']['mem_mb'],
        qos = res['run_DAS_Tool']['qos']
    benchmark:
        "output/benchmarks/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}_benchmark.txt"
    log:
        "output/logs/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}.log"
    shell:
        """
            DAS_Tool \
            --bins {input.metabat2},{input.maxbin2},{input.concoct} \
            --contigs {input.contigs} \
            --outputbasename {params.basename} \
            --labels metabat2,maxbin2,concoct \
            --write_bins \
            --debug \
            --write_bin_evals \
            --threads {threads} \
            --search_engine {params.search_engine} 2> {log} 1>&2
        """


rule consolidate_DAS_Tool_bins:
    """
    Consolidates and renames bin fastas generated by DAS_Tool into a single folder
    """
    input:
        "output/selected_bins/{mapper}/run_DAS_Tool/{contig_sample}/{contig_sample}_DASTool_summary.tsv"
    output:
        done = touch("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done")
    log:
        "output/logs/selected_bins/{mapper}/consolidate_DAS_Tool_bins/{contig_sample}.log"
    run:
        sample = wildcards.contig_sample 
        fasta_dir = join(dirname(input[0]),
                         sample + '_DASTool_bins')
        
        output_dir = dirname(output.done)
        # print(output_dir)
        # print(fasta_dir)
        fasta_files = glob(join(fasta_dir, '*.fa'))
        # print(fasta_files)
        for file in fasta_files:
            copyfile(file,
                     join(output_dir,
                          sample + '_' + basename(file)))

rule consolidate_DAS_Tool_bins_all:
    input:
        lambda wildcards: expand("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done",
                                 mapper=config['mappers'],
                                 contig_sample=contig_pairings.keys())


rule prepare_dRep:
    """
    Create file of paths for dRep
    """
    input:
        lambda wildcards: expand("output/selected_bins/{mapper}/DAS_Tool_Fastas/{contig_sample}.done",
                                 mapper=config['mappers'],
                                 contig_sample=contig_pairings.keys())
    output:
        "output/selected_bins/{mapper}/DAS_Tool_Fastas.input.txt"
    run:
        with open(output[0], 'w') as f:
            for p in input:
                f.write('%s\n' % p)


rule run_dRep:
    """
    Dereplicate bins using dRep
    """
    input:
         "output/selected_bins/{mapper}/DAS_Tool_Fastas.input.txt"
    output:
        outdir=directory("output/selected_bins/{mapper}/dRep"),
        outfig="output/selected_bins/{mapper}/dRep/figures/Winning_genomes.pdf"
    params:
        extra=config['params']['drep']['extra']
    conda:
        "../env/drep.yaml"
    threads:
        res['run_drep']['threads']
    resources:
        partition = res['run_drep']['partition'],
        mem_mb = res['run_drep']['mem_mb'],
        qos = res['run_drep']['qos']
    benchmark:
        "output/benchmarks/selected_bins/{mapper}/dRep/run_drep_benchmark.txt"
    log:
        "output/logs/selected_bins/{mapper}/dRep/run_drep.log"
    shell:
        """
            dRep dereplicate {output.outdir} {params.extra} \
              -p {threads} \
              -g {input} 2> {log} 1>&2
        """
