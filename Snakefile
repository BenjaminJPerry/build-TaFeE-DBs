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
        wget -O GTDB/bac120_metadata_latest.tar.gz {params.bacMeta};
        tar -xf GTDB/bac120_metadata_latest.tar.gz -O > {output.bac120Metadata};
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
        wget -O GTDB/ar53_metadata_latest.tar.gz {params.arcMeta};
        tar -xf GTDB/ar53_metadata_latest.tar.gz -O > {output.arc53Metadata};
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
        wget -O GTDB/bac120_taxonomy_latest.tsv.gz {params.gtdbBacTax};
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
        wget -O GTDB/ar53_taxonomy_latest.tsv.gz {params.gtdbArcTax};
        gunzip -c GTDB/ar53_taxonomy_latest.tsv.gz > {output.arcTax};
        '''


rule make_merged_taxonomy:
    output:
        metadata='GTDB/merged_taxonomy.tsv'
    input:
        host_taxonomy='resources/host_taxonomy.tsv',
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
        gtd_genomes_gz='GTDB/gtdb_genomes_reps_latest.tar.gz',
        sheep_gz = 'GTDB/host_genomes/sheep.fna.gz',
        cow_gz = 'GTDB/host_genomes/cow.fna.gz',
        goat_gz = 'GTDB/host_genomes/goat.fna.gz',
        deer_gz = 'GTDB/host_genomes/deer.fna.gz',
        wapiti_gz = 'GTDB/host_genomes/wapiti.fna.gz',
    threads: 2
    resources:
        partition='compute',
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
        mkdir GTDB/host_genomes
        wget -O {output.gtd_genomes_gz} {params.gtdbGenomes};
        wget -O {output.sheep_gz} {params.sheep};
        wget -O {output.cow_gz} {params.cow};
        wget -O {output.goat_gz} {params.goat};
        wget -O {output.deer_gz} {params.deer};
        wget -O {output.wapiti_gz} {params.wapiti};

        '''


rule prepare_GTDB_genomes:
    input:
        'GTDB/gtdb_genomes_reps_latest.tar.gz'
    output:
        directory('GTDB/input_genomes'),
    threads: 2
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 2 * 24 * 60
    shell:
        '''
        tar -xvzf {input};
        mkdir {output}; 
        find GTDB/gtdb_genomes_reps_latest -name "*.fna.gz" -exec mv -t {output}/ {} +;
        find GTDB/host_genomes -name "*.fna.gz" -exec mv -t {output}/ {} +;

        '''


rule prepKraken2Build:
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
        genomes_prep = directory('GTDB/kraken2-GTDB-214.1'),
        names_prep = 'GTDB/kraken2-GTDB-214.1/taxonomy/names.dmp',
        nodes_prep = 'GTDB/kraken2-GTDB-214.1/taxonomy/nodes.dmp',
    conda:
        'kraken2'
    threads: 64
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 5 * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 12
    shell:
        '''
        for file in $(ls {input.genomes});
        do
            kraken2-build --add-to-library {input.genomes}/$file --db {output}
        done

        mkdir -p {output}/taxonomy
        mv {input.nodes} {output}/taxonomy/{input.nodes}
        mv {input.names} {output}/taxonomy/{input.names}

        '''


rule build_kraken2:
    input:
        genomes_prep = 'GTDB/kraken2-GTDB-214.1',
        names_prep = 'GTDB/kraken2-GTDB-214.1/taxonomy/names.dmp',
        nodes_prep = 'GTDB/kraken2-GTDB-214.1/taxonomy/nodes.dmp',
    output:
        kraken2_index = 'GTDB/kraken2-GTDB-214.1/hash.k2d',
    conda:
        'kraken2'
    threads: 64
    resources:
        time = lambda wildcards, attempt: attempt * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 1600
    shell:
        '''
        kraken2-build --build --threads 64 --db {input}

        '''

#rule humann3_protein: #TODO

#rule humann3_default: #TODO
