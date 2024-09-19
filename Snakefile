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
        'GTDB214/kraken2-GTDB-214/hash.k2d',
        'GTDB214/merged_metadata.tsv',
        #'K2NT-20230205/hash.k2d',
        'GTDB214/kraken2-hosts/hash.k2d',


### Prepare ###
rule get_GTDB_bac_metadata:
    output:
        bac120Metadata='GTDB214/bac120_metadata_latest.tsv'
    threads: 2
    params:
        bacMeta=config['gtdb-bac-metadata']
    resources:
        partition='compute'
    shell:
        '''
        mkdir -p GTDB214
        curl -o GTDB214/bac120_metadata_latest.tsv.gz {params.bacMeta};
        gunzip -c GTDB214/bac120_metadata_latest.tsv.gz > {output.bac120Metadata};
        '''


rule get_GTDB_arc_metadata:
    output:
        arc53Metadata='GTDB214/ar53_metadata_latest.tsv'
    threads: 2
    params:
        arcMeta=config['gtdb-arc-metadata']
    resources:
        partition='compute'
    shell:
        '''
        mkdir -p GTDB214
        curl -o GTDB214/ar53_metadata_latest.tsv.gz {params.arcMeta};
        gunzip -c GTDB214/ar53_metadata_latest.tsv.gz > {output.arc53Metadata};
        '''


rule make_merged_metadata:
    output:
        metadata='GTDB214/merged_metadata.tsv'
    input:
        bac120Metadata='GTDB214/bac120_metadata_latest.tsv',
        arc53Metadata='GTDB214/ar53_metadata_latest.tsv',
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
        bacTax='GTDB214/bac120_taxonomy_latest.tsv'
    threads: 2
    resources:
        partition='compute'
    params:
        gtdbBacTax=config['gtdb-bac-tax']
    shell:
        '''
        mkdir -p GTDB214
        curl -o GTDB214/bac120_taxonomy_latest.tsv.gz {params.gtdbBacTax};
        gunzip -c GTDB214/bac120_taxonomy_latest.tsv.gz > {output.bacTax};
        '''


rule get_GTDB_arc_tax:
    output:
        arcTax='GTDB214/ar53_taxonomy_latest.tsv'
    threads: 2
    resources:
        partition='compute'
    params:
        gtdbArcTax = config['gtdb-arc-tax']
    shell:
        '''
        mkdir -p GTDB214
        curl -o GTDB214/ar53_taxonomy_latest.tsv.gz {params.gtdbArcTax};
        gunzip -c GTDB214/ar53_taxonomy_latest.tsv.gz > {output.arcTax};
        '''


rule make_merged_taxonomy:
    output:
        metadata='GTDB214/merged_taxonomy.tsv'
    input:
        host_taxonomy='resources/eukaryotic_taxa.tsv',
        arc_taxonomy='GTDB214/ar53_taxonomy_latest.tsv',
        bac_taxonomy='GTDB214/bac120_taxonomy_latest.tsv',
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
        gtdb_genomes_gz= protected('GTDB214/gtdb_genomes_reps_latest.tar.gz'),
        sheep_gz = protected('GTDB214/host_genomes/GCF_016772045.2_genomic.fna.gz'),
        cow_gz = protected('GTDB214/host_genomes/GCF_002263795.3_genomic.fna.gz'),
        goat_gz = protected('GTDB214/host_genomes/GCF_001704415.2_genomic.fna.gz'),
        deer_gz = protected('GTDB214/host_genomes/GCF_910594005.1_genomic.fna.gz'),
        wapiti_gz = protected('GTDB214/host_genomes/GCF_019320065.1_genomic.fna.gz'),
    threads: 2
    resources:
        partition='compute'
    resources:
        time = lambda wildcards, attempt: attempt * 1 * 24 * 60
    params:
        gtdbGenomes=config['gtdb-genomes'],
        sheep=config['sheep-genome'],
        cow=config['cow-genome'],
        goat=config['goat-genome'],
        deer=config['deer-genome'],
        wapiti=config['wapiti-genome'],
    shell:
        '''
        mkdir -p GTDB214/host_genomes

        curl -o {output.sheep_gz} {params.sheep};
        curl -o {output.cow_gz} {params.cow};
        curl -o {output.goat_gz} {params.goat};
        curl -o {output.deer_gz} {params.deer};
        curl -o {output.wapiti_gz} {params.wapiti};

        curl -o {output.gtdb_genomes_gz} {params.gtdbGenomes};

        '''


rule prepare_GTDB_genomes:
    input:
        gtdb= 'GTDB214/gtdb_genomes_reps_latest.tar.gz',
        sheep_gz = 'GTDB214/host_genomes/GCF_016772045.2_genomic.fna.gz',
        cow_gz = 'GTDB214/host_genomes/GCF_002263795.3_genomic.fna.gz',
        goat_gz = 'GTDB214/host_genomes/GCF_001704415.2_genomic.fna.gz',
        deer_gz = 'GTDB214/host_genomes/GCF_910594005.1_genomic.fna.gz',
        wapiti_gz = 'GTDB214/host_genomes/GCF_019320065.1_genomic.fna.gz',
    output:
        directory('GTDB214/input_genomes')
    threads: 2
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 1 * 24 * 60
    shell:
        '''

        mkdir -p {output}
        find GTDB214/host_genomes -name "*.fna.gz" -exec mv -t {output}/ {{}} +;

        mkdir GTDB214/gtdb_genomes_reps_latest
        tar -xvzf {input.gtdb} -C GTDB214/gtdb_genomes_reps_latest
        find GTDB214/gtdb_genomes_reps_latest -name "*.fna.gz" -exec mv -t {output}/ {{}} +;

        '''


rule prepare_kraken2_genomes:
    input:
        genomes='GTDB214/input_genomes',
        taxonomy='GTDB214/merged_taxonomy.tsv',
        tax_from_gtdb='workflow/scripts/tax_from_gtdb.py'
    output:
        genomes_out = directory('GTDB214/kraken_genomes'),
        nodes = 'GTDB214/nodes.dmp',
        names = 'GTDB214/names.dmp'
    conda:
        'kraken2'
    threads: 2
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 3 * 24 * 60,
    shell:
        '''
        python {input.tax_from_gtdb} --gtdb {input.taxonomy} --assemblies {input.genomes} --nodes {output.nodes} --names {output.names} --kraken_dir {output.genomes_out} &&
        rm -r GTDB214/gtdb_genomes_reps_latest
        rm -r GTDB214/host_genomes
        '''


rule prepare_kraken2_build:
    input:
        genomes = 'GTDB214/kraken_genomes',
        nodes = 'GTDB214/nodes.dmp',
        names = 'GTDB214/names.dmp',
    output:
        names_prep = 'GTDB214/kraken2-GTDB-214/taxonomy/names.dmp',
        nodes_prep = 'GTDB214/kraken2-GTDB-214/taxonomy/nodes.dmp',
    conda:
        'kraken2'
    threads: 32
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 5 * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 100
    shell:
        '''
        mkdir -p GTDB214/kraken2-GTDB-214/taxonomy

        cp {input.nodes} GTDB214/kraken2-GTDB-214/taxonomy/nodes.dmp
        cp {input.names} GTDB214/kraken2-GTDB-214/taxonomy/names.dmp

        for file in $(ls {input.genomes});
        do
            kraken2-build --threads {threads} --add-to-library {input.genomes}/$file --db GTDB214/kraken2-GTDB-214
        done


        '''


rule build_kraken2:
    input:
        names_prep = 'GTDB214/kraken2-GTDB-214/taxonomy/names.dmp',
        nodes_prep = 'GTDB214/kraken2-GTDB-214/taxonomy/nodes.dmp',
    output:
        kraken2_index = 'GTDB214/kraken2-GTDB-214/hash.k2d',
    conda:
        'kraken2'
    threads: 64
    benchmark:
        'benchmark/build_kraken2.txt'
    resources:
        partition='hugemem',
        time = lambda wildcards, attempt: attempt * 5 * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 1600
    shell:
        '''
        kraken2-build --build --threads {threads} --db GTDB214/kraken2-GTDB-214

        '''


rule prepare_kraken2_host_genomes:
    input:
        host_taxonomy='resources/host_taxa.tsv',
        tax_from_gtdb='workflow/scripts/tax_from_gtdb.py'
    output:
        genomes_out = directory('GTDB214/kraken_host_genomes'),
        nodes = 'GTDB214/nodes.host.dmp',
        names = 'GTDB214/names.host.dmp',
    conda:
        'kraken2'
    threads: 2
    resources:
        partition='compute',
        mem_gb = lambda wildcards, attempt: attempt * 16,
        time = lambda wildcards, attempt: attempt * 2 * 24 * 60,
    params:
        gtdbGenomes=config['gtdb-genomes'],
        sheep=config['sheep-genome'],
        cow=config['cow-genome'],
        goat=config['goat-genome'],
        deer=config['deer-genome'],
        wapiti=config['wapiti-genome'],
    shell:
        '''
        mkdir -p GTDB214/host_genomes

        curl -o GTDB214/host_genomes/GCF_016772045.2_genomic.fna.gz {params.sheep};
        curl -o GTDB214/host_genomes/GCF_002263795.3_genomic.fna.gz {params.cow};
        curl -o GTDB214/host_genomes/GCF_001704415.2_genomic.fna.gz {params.goat};
        curl -o GTDB214/host_genomes/GCF_910594005.1_genomic.fna.gz {params.deer};
        curl -o GTDB214/host_genomes/GCF_019320065.1_genomic.fna.gz {params.wapiti};


        python {input.tax_from_gtdb} --gtdb {input.host_taxonomy} --assemblies GTDB214/host_genomes --nodes {output.nodes} --names {output.names} --kraken_dir {output.genomes_out}

        '''


rule prepare_kraken2_hosts_build:
    input:
        genomes = 'GTDB214/kraken_host_genomes',
        nodes = 'GTDB214/nodes.host.dmp',
        names = 'GTDB214/names.host.dmp'
    output:
        names_prep = 'GTDB214/kraken2-hosts/taxonomy/names.dmp',
        nodes_prep = 'GTDB214/kraken2-hosts/taxonomy/nodes.dmp',
    conda:
        'kraken2'
    threads: 16
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 2 * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 24
    shell:
        '''
        mkdir -p GTDB214/kraken2-hosts/taxonomy

        cp {input.nodes} GTDB214/kraken2-hosts/taxonomy/nodes.dmp
        cp {input.names} GTDB214/kraken2-hosts/taxonomy/names.dmp

        for file in $(ls {input.genomes});
        do
            kraken2-build --threads {threads} --add-to-library {input.genomes}/$file --db GTDB214/kraken2-hosts
        done


        '''


rule build_kraken2_hosts:
    input:
        names_prep = 'GTDB214/kraken2-hosts/taxonomy/names.dmp',
        nodes_prep = 'GTDB214/kraken2-hosts/taxonomy/nodes.dmp',
    output:
        kraken2_index = 'GTDB214/kraken2-hosts/hash.k2d',
    conda:
        'kraken2'
    threads: 64
    benchmark:
        'benchmark/build_kraken2_hosts.txt'
    resources:
        partition='hugemem',
        time = lambda wildcards, attempt: attempt * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 1600
    shell:
        '''
        kraken2-build --build --threads {threads} --db GTDB214/kraken2-hosts
        '''


