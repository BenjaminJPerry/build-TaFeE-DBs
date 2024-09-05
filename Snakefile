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
        'GTDB/kraken2-GTDB-220.0/hash.k2d',
        'GTDB/merged_metadata.tsv',
        #'K2NT-20230205/hash.k2d',
        'GTDB/kraken2-hosts/hash.k2d',


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
        curl -o GTDB/bac120_metadata_latest.tsv.gz {params.bacMeta};
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
        curl -o GTDB/ar53_metadata_latest.tsv.gz {params.arcMeta};
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
        curl -o GTDB/bac120_taxonomy_latest.tsv.gz {params.gtdbBacTax};
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
        curl -o GTDB/ar53_taxonomy_latest.tsv.gz {params.gtdbArcTax};
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
        gtdb_genomes_gz= protected('GTDB/gtdb_genomes_reps_latest.tar.gz'),
        sheep_gz = protected('GTDB/host_genomes/GCF_016772045.2_genomic.fna.gz'),
        cow_gz = protected('GTDB/host_genomes/GCF_002263795.3_genomic.fna.gz'),
        goat_gz = protected('GTDB/host_genomes/GCF_001704415.2_genomic.fna.gz'),
        deer_gz = protected('GTDB/host_genomes/GCF_910594005.1_genomic.fna.gz'),
        wapiti_gz = protected('GTDB/host_genomes/GCF_019320065.1_genomic.fna.gz'),
        Entodinium_caudatum_gz = protected('GTDB/protist_genomes/GCA_002087855.3_genomic.fna.gz'),
        Entodinium_longinucleatum_gz = protected('GTDB/protist_genomes/GCA_023897235.1_genomic.fna.gz'),
        Entodinium_bursa_gz = protected('GTDB/protist_genomes/GCA_023807345.1_genomic.fna.gz'),
        Epidinium_caudatum_gz = protected('GTDB/protist_genomes/GCA_023807225.1_genomic.fna.gz'),
        Epidinium_cattanei_gz = protected('GTDB/protist_genomes/GCA_023805625.1_genomic.fna.gz'),
        Diplodinium_dentatum_gz = protected('GTDB/protist_genomes/GCA_023807165.1_genomic.fna.gz'),
        Diplodinium_flabellum_gz = protected('GTDB/protist_genomes/GCA_023806845.1_genomic.fna.gz'),
        Isotricha_intestinalis_gz = protected('GTDB/protist_genomes/GCA_023807065.1_genomic.fna.gz'),
        Isotricha_prostoma_gz = protected('GTDB/protist_genomes/GCA_023807205.1_genomic.fna.gz'),
        Isotricha_YL_2021a_gz = protected('GTDB/protist_genomes/GCA_023806865.1_genomic.fna.gz'),
        Isotricha_YL_2021b_gz = protected('GTDB/protist_genomes/GCA_023805745.1_genomic.fna.gz'),
        Dasytricha_ruminantium_gz = protected('GTDB/protist_genomes/GCA_023805585.1_genomic.fna.gz'),
        Ophryoscolex_caudatus_gz = protected('GTDB/protist_genomes/GCA_023806825.1_genomic.fna.gz'),
        Polyplastron_multivesiculatum_gz= protected('GTDB/protist_genomes/GCA_023783355.1_genomic.fna.gz'),
        Eremoplastron_rostratum_gz = protected('GTDB/protist_genomes/GCA_023805755.1_genomic.fna.gz'),
        Ostracodinium_gracile_gz = protected('GTDB/protist_genomes/GCA_023805685.1_genomic.fna.gz'),
        Metadinium_minorum_gz = protected('GTDB/protist_genomes/GCA_023807265.1_genomic.fna.gz'),
        Enoploplastron_triloricatum_gz = protected('GTDB/protist_genomes/GCA_023783335.1_genomic.fna.gz'),
        Ostracodinium_dentatum_gz = protected('GTDB/protist_genomes/GCA_023805525.1_genomic.fna.gz'),
        Anaeromyces_S4_gz = protected('GTDB/fungi_genomes/GCA_002104895.1_genomic.fna.gz'),
        Caecomyces_SIG737_gz = protected('GTDB/fungi_genomes/GCA_027248405.1_genomic.fna.gz'),
        Neocallimastix_californiae_gz = protected('GTDB/fungi_genomes/GCA_002104975.1_genomic.fna.gz'),
        Neocallimastix_JGI_2020a_gz = protected('GTDB/fungi_genomes/GCA_016946835.1_genomic.fna.gz'),
        Piromyces_finnis_gz = protected('GTDB/fungi_genomes/GCA_002104945.1_genomic.fna.gz'),
        Piromyces_E2_gz = protected('GTDB/fungi_genomes/GCA_002157105.1_genomic.fna.gz'),
        Piromyces_SIG733_gz = protected('GTDB/fungi_genomes/GCA_027248865.1_genomic.fna.gz'),
        Orpinomyces_sp_gz = protected('GTDB/fungi_genomes/GCA_000412615.1_genomic.fna.gz'),
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
        Entodinium_caudatum=config['Entodinium-caudatum'],
        Entodinium_longinucleatum=config['Entodinium-longinucleatum'],
        Entodinium_bursa=config['Entodinium-bursa'],
        Epidinium_caudatum=config['Epidinium-caudatum'],
        Epidinium_cattanei=config['Epidinium-cattanei'],
        Diplodinium_dentatum=config['Diplodinium-dentatum'],
        Diplodinium_flabellum=config['Diplodinium-flabellum'],
        Isotricha_intestinalis=config['Isotricha-intestinalis'],
        Isotricha_prostoma=config['Isotricha-prostoma'],
        Isotricha_YL_2021a=config['Isotricha-YL-2021a'],
        Isotricha_YL_2021b=config['Isotricha-YL-2021b'],
        Dasytricha_ruminantium=config['Dasytricha-ruminantium'],
        Ophryoscolex_caudatus=config['Ophryoscolex-caudatus'],
        Polyplastron_multivesiculatum=config['Polyplastron-multivesiculatum'],
        Eremoplastron_rostratum=config['Eremoplastron-rostratum'],
        Ostracodinium_gracile=config['Ostracodinium-gracile'],
        Metadinium_minorum=config['Metadinium-minorum'],
        Enoploplastron_triloricatum=config['Enoploplastron-triloricatum'],
        Ostracodinium_dentatum=config['Ostracodinium-dentatum'],
        Anaeromyces_S4=config['Anaeromyces-S4'],
        Caecomyces_SIG737=config['Caecomyces-SIG737'],
        Neocallimastix_californiae=config['Neocallimastix-californiae'],
        Neocallimastix_JGI_2020a=config['Neocallimastix-JGI-2020a'],
        Piromyces_finnis=config['Piromyces-finnis'],
        Piromyces_E2=config['Piromyces-E2'],
        Piromyces_SIG733=config['Piromyces-SIG733'],
        Orpinomyces_sp=config['Orpinomyces-sp'],
    shell:
        '''
        mkdir -p GTDB/host_genomes
        mkdir -p GTDB/protist_genomes
        mkdir -p GTDB/fungi_genomes

        curl -o {output.sheep_gz} {params.sheep};
        curl -o {output.cow_gz} {params.cow};
        curl -o {output.goat_gz} {params.goat};
        curl -o {output.deer_gz} {params.deer};
        curl -o {output.wapiti_gz} {params.wapiti};

        curl -o {output.Entodinium_caudatum_gz} {params.Entodinium_caudatum};
        curl -o {output.Entodinium_longinucleatum_gz} {params.Entodinium_longinucleatum};
        curl -o {output.Entodinium_bursa_gz} {params.Entodinium_bursa};
        curl -o {output.Epidinium_caudatum_gz} {params.Epidinium_caudatum};
        curl -o {output.Epidinium_cattanei_gz} {params.Epidinium_cattanei};
        curl -o {output.Diplodinium_dentatum_gz} {params.Diplodinium_dentatum};
        curl -o {output.Diplodinium_flabellum_gz} {params.Diplodinium_flabellum};
        curl -o {output.Isotricha_intestinalis_gz} {params.Isotricha_intestinalis};
        curl -o {output.Isotricha_prostoma_gz} {params.Isotricha_prostoma};
        curl -o {output.Isotricha_YL_2021a_gz} {params.Isotricha_YL_2021a};
        curl -o {output.Isotricha_YL_2021b_gz} {params.Isotricha_YL_2021b};
        curl -o {output.Dasytricha_ruminantium_gz} {params.Dasytricha_ruminantium};
        curl -o {output.Ophryoscolex_caudatus_gz} {params.Ophryoscolex_caudatus};
        curl -o {output.Polyplastron_multivesiculatum_gz} {params.Polyplastron_multivesiculatum};
        curl -o {output.Eremoplastron_rostratum_gz} {params.Eremoplastron_rostratum};
        curl -o {output.Ostracodinium_gracile_gz} {params.Ostracodinium_gracile};
        curl -o {output.Metadinium_minorum_gz} {params.Metadinium_minorum};
        curl -o {output.Enoploplastron_triloricatum_gz} {params.Enoploplastron_triloricatum};
        curl -o {output.Ostracodinium_dentatum_gz} {params.Ostracodinium_dentatum};
        curl -o {output.Anaeromyces_S4_gz} {params.Anaeromyces_S4};
        curl -o {output.Caecomyces_SIG737_gz} {params.Caecomyces_SIG737};
        curl -o {output.Neocallimastix_californiae_gz} {params.Neocallimastix_californiae};
        curl -o {output.Neocallimastix_JGI_2020a_gz} {params.Neocallimastix_JGI_2020a};
        curl -o {output.Piromyces_finnis_gz} {params.Piromyces_finnis};
        curl -o {output.Piromyces_E2_gz} {params.Piromyces_E2};
        curl -o {output.Piromyces_SIG733_gz} {params.Piromyces_SIG733};
        curl -o {output.Orpinomyces_sp_gz} {params.Orpinomyces_sp};

        curl -o {output.gtdb_genomes_gz} {params.gtdbGenomes};

        '''


rule prepare_GTDB_genomes:
    input:
        gtdb= 'GTDB/gtdb_genomes_reps_latest.tar.gz',
        sheep_gz = 'GTDB/host_genomes/GCF_016772045.2_genomic.fna.gz',
        cow_gz = 'GTDB/host_genomes/GCF_002263795.3_genomic.fna.gz',
        goat_gz = 'GTDB/host_genomes/GCF_001704415.2_genomic.fna.gz',
        deer_gz = 'GTDB/host_genomes/GCF_910594005.1_genomic.fna.gz',
        wapiti_gz = 'GTDB/host_genomes/GCF_019320065.1_genomic.fna.gz',
        Entodinium_caudatum_gz = 'GTDB/protist_genomes/GCA_002087855.3_genomic.fna.gz',
        Entodinium_longinucleatum_gz = 'GTDB/protist_genomes/GCA_023897235.1_genomic.fna.gz',
        Entodinium_bursa_gz = 'GTDB/protist_genomes/GCA_023807345.1_genomic.fna.gz',
        Epidinium_caudatum_gz = 'GTDB/protist_genomes/GCA_023807225.1_genomic.fna.gz',
        Epidinium_cattanei_gz = 'GTDB/protist_genomes/GCA_023805625.1_genomic.fna.gz',
        Diplodinium_dentatum_gz = 'GTDB/protist_genomes/GCA_023807165.1_genomic.fna.gz',
        Diplodinium_flabellum_gz = 'GTDB/protist_genomes/GCA_023806845.1_genomic.fna.gz',
        Isotricha_intestinalis_gz = 'GTDB/protist_genomes/GCA_023807065.1_genomic.fna.gz',
        Isotricha_prostoma_gz = 'GTDB/protist_genomes/GCA_023807205.1_genomic.fna.gz',
        Isotricha_YL_2021a_gz = 'GTDB/protist_genomes/GCA_023806865.1_genomic.fna.gz',
        Isotricha_YL_2021b_gz = 'GTDB/protist_genomes/GCA_023805745.1_genomic.fna.gz',
        Dasytricha_ruminantium_gz = 'GTDB/protist_genomes/GCA_023805585.1_genomic.fna.gz',
        Ophryoscolex_caudatus_gz = 'GTDB/protist_genomes/GCA_023806825.1_genomic.fna.gz',
        Polyplastron_multivesiculatum_gz= 'GTDB/protist_genomes/GCA_023783355.1_genomic.fna.gz',
        Eremoplastron_rostratum_gz = 'GTDB/protist_genomes/GCA_023805755.1_genomic.fna.gz',
        Ostracodinium_gracile_gz = 'GTDB/protist_genomes/GCA_023805525.1_genomic.fna.gz',
        Metadinium_minorum_gz = 'GTDB/protist_genomes/GCA_023807265.1_genomic.fna.gz',
        Enoploplastron_triloricatum_gz = 'GTDB/protist_genomes/GCA_023783335.1_genomic.fna.gz',
        Ostracodinium_dentatum_gz = 'GTDB/protist_genomes/GCA_023805525.1_genomic.fna.gz',
        Anaeromyces_S4_gz = 'GTDB/fungi_genomes/GCA_002104895.1_genomic.fna.gz',
        Caecomyces_SIG737_gz = 'GTDB/fungi_genomes/GCA_027248405.1_genomic.fna.gz',
        Neocallimastix_californiae_gz = 'GTDB/fungi_genomes/GCA_002104975.1_genomic.fna.gz',
        Neocallimastix_JGI_2020a_gz = 'GTDB/fungi_genomes/GCA_016946835.1_genomic.fna.gz',
        Piromyces_finnis_gz = 'GTDB/fungi_genomes/GCA_002104945.1_genomic.fna.gz',
        Piromyces_E2_gz = 'GTDB/fungi_genomes/GCA_002157105.1_genomic.fna.gz',
        Piromyces_SIG733_gz = 'GTDB/fungi_genomes/GCA_027248865.1_genomic.fna.gz',
        Orpinomyces_sp_gz = 'GTDB/fungi_genomes/GCA_000412615.1_genomic.fna.gz',
    output:
        directory('GTDB/input_genomes')
    threads: 2
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 1 * 24 * 60
    shell:
        '''

        mkdir -p {output}
        find GTDB/host_genomes -name "*.fna.gz" -exec mv -t {output}/ {{}} +;

        find GTDB/protist_genomes -name "*.fna.gz" -exec mv -t {output}/ {{}} +;
        find GTDB/fungi_genomes -name "*.fna.gz" -exec mv -t {output}/ {{}} +;

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
        time = lambda wildcards, attempt: attempt * 3 * 24 * 60,
    shell:
        '''
        python {input.tax_from_gtdb} --gtdb {input.taxonomy} --assemblies {input.genomes} --nodes {output.nodes} --names {output.names} --kraken_dir {output.genomes_out} &&
        rm -r GTDB/gtdb_genomes_reps_latest
        rm -r GTDB/host_genomes
        rm -r GTDB/protist_genomes
        rm -r GTDB/fungi_genomes
        '''


rule prepare_kraken2_build:
    input:
        genomes = 'GTDB/kraken_genomes',
        nodes = 'GTDB/nodes.dmp',
        names = 'GTDB/names.dmp',
    output:
        names_prep = 'GTDB/kraken2-GTDB-220.0/taxonomy/names.dmp',
        nodes_prep = 'GTDB/kraken2-GTDB-220.0/taxonomy/nodes.dmp',
    conda:
        'kraken2'
    threads: 32
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 5 * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 100
    shell:
        '''
        mkdir -p GTDB/kraken2-GTDB-220.0/taxonomy

        cp {input.nodes} GTDB/kraken2-GTDB-220.0/taxonomy/nodes.dmp
        cp {input.names} GTDB/kraken2-GTDB-220.0/taxonomy/names.dmp

        for file in $(ls {input.genomes});
        do
            kraken2-build --threads {threads} --add-to-library {input.genomes}/$file --db GTDB/kraken2-GTDB-220.0
        done


        '''


rule build_kraken2:
    input:
        names_prep = 'GTDB/kraken2-GTDB-220.0/taxonomy/names.dmp',
        nodes_prep = 'GTDB/kraken2-GTDB-220.0/taxonomy/nodes.dmp',
    output:
        kraken2_index = 'GTDB/kraken2-GTDB-220.0/hash.k2d',
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
        kraken2-build --build --threads {threads} --db GTDB/kraken2-GTDB-220.0

        '''


rule prepare_kraken2_host_genomes:
    input:
        host_taxonomy='resources/eukaryotic_taxa.tsv',
        tax_from_gtdb='workflow/scripts/tax_from_gtdb.py'
    output:
        genomes_out = directory('GTDB/kraken_host_genomes'),
        nodes = 'GTDB/nodes.host.dmp',
        names = 'GTDB/names.host.dmp',
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
        mkdir -p GTDB/host_genomes

        curl -o GTDB/host_genomes/GCF_016772045.2_genomic.fna.gz {params.sheep};
        curl -o GTDB/host_genomes/GCF_002263795.3_genomic.fna.gz {params.cow};
        curl -o GTDB/host_genomes/GCF_001704415.2_genomic.fna.gz {params.goat};
        curl -o GTDB/host_genomes/GCF_910594005.1_genomic.fna.gz {params.deer};
        curl -o GTDB/host_genomes/GCF_019320065.1_genomic.fna.gz {params.wapiti};


        python {input.tax_from_gtdb} --gtdb {input.host_taxonomy} --assemblies GTDB/host_genomes --nodes {output.nodes} --names {output.names} --kraken_dir {output.genomes_out}

        '''


rule prepare_kraken2_hosts_build:
    input:
        genomes = 'GTDB/kraken_host_genomes',
        nodes = 'GTDB/nodes.host.dmp',
        names = 'GTDB/names.host.dmp'
    output:
        names_prep = 'GTDB/kraken2-hosts/taxonomy/names.dmp',
        nodes_prep = 'GTDB/kraken2-hosts/taxonomy/nodes.dmp',
    conda:
        'kraken2'
    threads: 16
    resources:
        partition='compute',
        time = lambda wildcards, attempt: attempt * 2 * 24 * 60,
        mem_gb = lambda wildcards, attempt: attempt * 24
    shell:
        '''
        mkdir -p GTDB/kraken2-hosts/taxonomy

        cp {input.nodes} GTDB/kraken2-hosts/taxonomy/nodes.dmp
        cp {input.names} GTDB/kraken2-hosts/taxonomy/names.dmp

        for file in $(ls {input.genomes});
        do
            kraken2-build --threads {threads} --add-to-library {input.genomes}/$file --db GTDB/kraken2-hosts
        done


        '''


rule build_kraken2_hosts:
    input:
        names_prep = 'GTDB/kraken2-hosts/taxonomy/names.dmp',
        nodes_prep = 'GTDB/kraken2-hosts/taxonomy/nodes.dmp',
    output:
        kraken2_index = 'GTDB/kraken2-hosts/hash.k2d',
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
        kraken2-build --build --threads {threads} --db GTDB/kraken2-hosts
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
        curl -o K2NT-20230205.tar.gz {params.k2_prebuilt_nt};
        mkdir -p K2NT-20230205
        tar -xvzf K2NT-20230205.tar.gz -C K2NT-20230205
        '''


#rule humann3_protein: #TODO
    

#rule humann3_default: #TODO
