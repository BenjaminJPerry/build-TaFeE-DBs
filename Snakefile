# Copyright (c) 2023 Benjamin J Perry
# Version: alpha
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"

import os

onstart:
    print(f"Working directory: {os.getcwd()}")

    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")

    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')


rule targets:
    input:
        'GTDB/kraken2-GTDB-214.1/hash.k2d',
        'GTDB/merged_metadata.tsv',
        'K2NT-20230205/hash.k2d',
        #'biobakery/...' #TODO


### Prepare ###
rule get_GTDB_bac_metadata:
    output:
        bac120Metadata='GTDB/bac120_metadata_latest.tsv'
    threads: 2
    params:
        bacMeta=config['gtdb-bac-metadata']
    resources:
        partition='compute'
    shell:
        '''
        mkdir -p GTDB
        wget -c -O GTDB/bac120_metadata_latest.tsv.gz {params.bacMeta};
        gunzip -c GTDB/bac120_metadata_latest.tsv.gz > {output.bac120Metadata};
        '''


rule get_GTDB_arc_metadata:
    output:
        arc53Metadata='GTDB/ar53_metadata_latest.tsv'
    threads: 2
    params:
        arcMeta=config['gtdb-arc-metadata']
    resources:
        partition='compute'
    shell:
        '''
        mkdir -p GTDB
        wget -c -O GTDB/ar53_metadata_latest.tsv.gz {params.arcMeta};
        gunzip -c GTDB/ar53_metadata_latest.tsv.gz > {output.arc53Metadata};
        '''


rule make_merged_metadata:
    output:
        metadata='GTDB/merged_metadata.tsv'
    input:
        bac120Metadata='GTDB/bac120_metadata_latest.tsv',
        arc53Metadata='GTDB/ar53_metadata_latest.tsv',
    threads: 2
    resources:
        partition='compute'
    shell:
        '''
        cat {input.bac120Metadata} > {output.metadata};
        cat {input.arc53Metadata} | grep -v "accession" >> {output.metadata};
        '''


rule get_GTDB_bac_tax:
    output:
        bacTax='GTDB/bac120_taxonomy_latest.tsv'
    threads: 2
    resources:
        partition='compute'
    params:
        gtdbBacTax=config['gtdb-bac-tax']
    shell:
        '''
        mkdir -p GTDB
        wget -c -O GTDB/bac120_taxonomy_latest.tsv.gz {params.gtdbBacTax};
        gunzip -c GTDB/bac120_taxonomy_latest.tsv.gz > {output.bacTax};
        '''


rule get_GTDB_arc_tax:
    output:
        arcTax='GTDB/ar53_taxonomy_latest.tsv'
    threads: 2
    resources:
        partition='compute'
    params:
        gtdbArcTax = config['gtdb-arc-tax']
    shell:
        '''
        mkdir -p GTDB
        wget -c -O GTDB/ar53_taxonomy_latest.tsv.gz {params.gtdbArcTax};
        gunzip -c GTDB/ar53_taxonomy_latest.tsv.gz > {output.arcTax};
        '''


rule make_merged_taxonomy:
    output:
        metadata='GTDB/merged_taxonomy.tsv'
    input:
        host_taxonomy='resources/eukaryotic_taxa.tsv',
        arc_taxonomy='GTDB/ar53_taxonomy_latest.tsv',
        bac_taxonomy='GTDB/bac120_taxonomy_latest.tsv',
    threads: 2
    resources:
        partition='compute'
    shell:
        '''
        cat {input.host_taxonomy} > {output.metadata};
        cat {input.arc_taxonomy} | grep -v "accession" >> {output.metadata};
        cat {input.bac_taxonomy} | grep -v "accession" >> {output.metadata};
        '''


rule get_genomes:
    output:
        gtd_genomes_gz= protected('GTDB/gtdb_genomes_reps_latest.tar.gz'),
        sheep_gz = protected('GTDB/host_genomes/GCF_000298735.2_genomic.fna.gz'),
        cow_gz = protected('GTDB/host_genomes/GCF_002263795.3_genomic.fna.gz'),
        goat_gz = protected('GTDB/host_genomes/GCF_001704415.2_genomic.fna.gz'),
        deer_gz = protected('GTDB/host_genomes/GCF_910594005.1_genomic.fna.gz'),
        wapiti_gz = protected('GTDB/host_genomes/GCF_019320065.1_genomic.fna.gz'),
    threads: 2
    resources:
        partition='compute'
    resources:
        time = lambda wildcards, attempt: attempt * 7 * 24 * 60
    params:
        gtdbGenomes=config['gtdb-genomes'],
        sheep=config['sheep-genome'],
        cow=config['cow-genome'],
        goat=config['goat-genome'],
        deer=config['deer-genome'],
        wapiti=config['wapiti-genome'],
    shell:
        '''
        mkdir -p GTDB/host_genomes

        wget -c -O {output.sheep_gz} {params.sheep};
        wget -c -O {output.cow_gz} {params.cow};
        wget -c -O {output.goat_gz} {params.goat};
        wget -c -O {output.deer_gz} {params.deer};
        wget -c -O {output.wapiti_gz} {params.wapiti};

        wget -c -O {output.gtd_genomes_gz} {params.gtdbGenomes};

        '''


rule prepare_GTDB_genomes:
    input:
        gtdb = 'GTDB/gtdb_genomes_reps_latest.tar.gz',
        sheep_gz = 'GTDB/host_genomes/GCF_000298735.2_genomic.fna.gz',
        cow_gz = 'GTDB/host_genomes/GCF_002263795.3_genomic.fna.gz',
        goat_gz = 'GTDB/host_genomes/GCF_001704415.2_genomic.fna.gz',
        deer_gz = 'GTDB/host_genomes/GCF_910594005.1_genomic.fna.gz',
        wapiti_gz = 'GTDB/host_genomes/GCF_019320065.1_genomic.fna.gz',
    output:
        directory('GTDB/input_genomes')
    threads: 2
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 2 * 24 * 60
    shell:
        '''

        mkdir -p {output}
        find GTDB/host_genomes -name "*.fna.gz" -exec mv -t {output}/ {{}} +;

        mkdir GTDB/gtdb_genomes_reps_latest
        tar -xvzf {input.gtdb} -C GTDB/gtdb_genomes_reps_latest
        find GTDB/gtdb_genomes_reps_latest -name "*.fna.gz" -exec mv -t {output}/ {{}} +;

        '''


rule prepare_kraken2_genomes:
    input:
        genomes='GTDB/input_genomes',
        taxonomy='GTDB/merged_taxonomy.tsv',
        tax_from_gtdb='workflow/scripts/tax_from_gtdb.py'
    output:
        genomes_out = directory('GTDB/kraken_genomes'),
        nodes = 'GTDB/nodes.dmp',
        names = 'GTDB/names.dmp'
    conda:
        'kraken2'
    threads: 2
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 2 * 24 * 60,
    shell:
        '''
        python {input.tax_from_gtdb} --gtdb {input.taxonomy} --assemblies {input.genomes} --nodes {output.nodes} --names {output.names} --kraken_dir {output.genomes_out} &&
        rm -r GTDB/gtdb_genomes_reps_latest
        rm -r GTDB/host_genomes

        '''


rule prepare_kraken2_build:
    input:
        genomes = 'GTDB/kraken_genomes',
        nodes = 'GTDB/nodes.dmp',
        names = 'GTDB/names.dmp',
    output:
        names_prep = 'GTDB/kraken2-GTDB-214.1/taxonomy/names.dmp',
        nodes_prep = 'GTDB/kraken2-GTDB-214.1/taxonomy/nodes.dmp',
    conda:
        'kraken2'
    threads: 16
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 5 * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 24
    shell:
        '''
        mkdir -p GTDB/kraken2-GTDB-214.1/taxonomy

        cp {input.nodes} GTDB/kraken2-GTDB-214.1/taxonomy/nodes.dmp
        cp {input.names} GTDB/kraken2-GTDB-214.1/taxonomy/names.dmp

        for file in $(ls {input.genomes});
        do
            kraken2-build --threads {threads} --add-to-library {input.genomes}/$file --db GTDB/kraken2-GTDB-214.1
        done


        '''


rule build_kraken2:
    input:
        names_prep = 'GTDB/kraken2-GTDB-214.1/taxonomy/names.dmp',
        nodes_prep = 'GTDB/kraken2-GTDB-214.1/taxonomy/nodes.dmp',
    output:
        kraken2_index = 'GTDB/kraken2-GTDB-214.1/hash.k2d',
    conda:
        'kraken2'
    threads: 64
    resources:
        partition='hugemem',
        time = lambda wildcards, attempt: attempt * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 1600
    shell:
        '''
        kraken2-build --build --threads {threads} --db GTDB/kraken2-GTDB-214.1

        '''


rule kraken2_prebuilt_ntdb:
    output:
        'K2NT-20230205/hash.k2d'
    threads: 2
    resources:
        partition='compute'
    params:
        k2_prebuilt_nt = config['k2nt']
    shell:
        '''
        wget -c -O K2NT-20230205.tar.gz {params.k2_prebuilt_nt};
        mkdir -p K2NT-20230205
        tar -xvzf K2NT-20230205.tar.gz -C K2NT-20230205
        '''


#rule humann3_protein: #TODO
    

#rule humann3_default: #TODO
