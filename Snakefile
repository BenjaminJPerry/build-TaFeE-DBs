# Copyright (c) 2023-2026 Benjamin J Perry
# Generalized kraken2 build pipeline: user supplies (accessions.txt, optional GTDB release)
# and the workflow auto-derives URLs, taxonomy, and a provenance-named output directory.
#
# Invoke with: snakemake --configfile config/config.yaml ...
# (No default configfile is loaded — pass it explicitly to avoid accidental merges.)

import os
import sys
from pathlib import Path

if not config:
    sys.exit("Error: no config loaded. Pass --configfile config/config.yaml (or config/test_config.yaml).")

sys.path.insert(0, "workflow/scripts")
from compute_db_name import read_accessions, canonical, compute_name
from resolve_gtdb_urls import gtdb_urls


onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS:")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')


OUTPUT_ROOT = config.get("output_root", "GTDB")
DATABASES = config.get("databases", {})
if not DATABASES:
    sys.exit("Error: config must define a 'databases' dict")

# DAG-time resolution: config-key -> (directory name, accessions list, gtdb URLs or None)
DB_DIRS = {}
DB_ACCESSIONS = {}
DB_BY_DIR = {}
DB_GTDB_URLS = {}
DB_SPEC = {}

for key, spec in DATABASES.items():
    acc_path = Path(spec["supplemental_accessions"])
    if not acc_path.exists():
        sys.exit(f"Error: accessions file missing for db '{key}': {acc_path}")
    accs = read_accessions(acc_path)
    release = spec.get("gtdb_release")
    label = spec.get("label", key)
    name = spec.get("name_override") or compute_name(label, accs, release)
    DB_DIRS[key] = name
    DB_ACCESSIONS[key] = accs
    DB_BY_DIR[name] = key
    DB_GTDB_URLS[key] = gtdb_urls(str(release)) if release else None
    DB_SPEC[key] = spec
    print(f"  db '{key}' -> {OUTPUT_ROOT}/{name}  ({len(accs)} accessions, gtdb={release or 'none'})")


def spec_of(db_dir: str):
    return DB_SPEC[DB_BY_DIR[db_dir]]


def accessions_of(db_dir: str):
    return DB_ACCESSIONS[DB_BY_DIR[db_dir]]


def has_gtdb(db_dir: str) -> bool:
    return DB_GTDB_URLS[DB_BY_DIR[db_dir]] is not None


def gtdb_url(db_dir: str, key: str) -> str:
    urls = DB_GTDB_URLS[DB_BY_DIR[db_dir]]
    if not urls:
        sys.exit(f"db '{db_dir}' has no gtdb_release configured")
    return urls[key]


def res(db_dir: str, rule_key: str, field: str, default):
    return spec_of(db_dir).get("resources", {}).get(rule_key, {}).get(field, default)


def merged_tax_inputs(wildcards):
    out = {"supp": f"{OUTPUT_ROOT}/{wildcards.db}/supplemental.taxonomy.tsv"}
    if has_gtdb(wildcards.db):
        out["bac"] = f"{OUTPUT_ROOT}/{wildcards.db}/gtdb/bac120_taxonomy.tsv"
        out["arc"] = f"{OUTPUT_ROOT}/{wildcards.db}/gtdb/ar53_taxonomy.tsv"
    return out


wildcard_constraints:
    db = r"kraken2-[A-Za-z0-9._-]+",
    accession = r"(?:GCA|GCF)_\d{9}\.\d+",


rule targets:
    input:
        expand(f"{OUTPUT_ROOT}/{{db}}/hash.k2d", db=DB_DIRS.values()),
        expand(f"{OUTPUT_ROOT}/{{db}}/MANIFEST.{{db}}.json", db=DB_DIRS.values()),


# ---- Supplemental taxonomy + URL resolution (always runs) -------------------

# Shared taxdump dir (NCBI nodes.dmp/names.dmp/merged.dmp/delnodes.dmp) used by taxonkit
# to resolve NCBI taxids into lineage strings. Downloaded once per output_root.
TAXDUMP_DIR = f"{OUTPUT_ROOT}/_taxdump"


rule download_taxdump:
    output:
        nodes=f"{TAXDUMP_DIR}/nodes.dmp",
        names=f"{TAXDUMP_DIR}/names.dmp",
        merged=f"{TAXDUMP_DIR}/merged.dmp",
        delnodes=f"{TAXDUMP_DIR}/delnodes.dmp",
    threads: 2
    resources: mem_gb=4, time=60, partition="compute"
    shell:
        """
        mkdir -p {TAXDUMP_DIR}
        curl -fsSL --retry 5 --retry-delay 15 \
            -o {TAXDUMP_DIR}/taxdump.tar.gz \
            https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -xzf {TAXDUMP_DIR}/taxdump.tar.gz -C {TAXDUMP_DIR} \
            nodes.dmp names.dmp merged.dmp delnodes.dmp
        rm -f {TAXDUMP_DIR}/taxdump.tar.gz
        """


rule fetch_supplemental_taxonomy:
    output:
        tax=f"{OUTPUT_ROOT}/{{db}}/supplemental.taxonomy.tsv",
        urls=f"{OUTPUT_ROOT}/{{db}}/supplemental.urls.tsv",
    input:
        taxdump_nodes=f"{TAXDUMP_DIR}/nodes.dmp",
    params:
        accessions=lambda wc: spec_of(wc.db)["supplemental_accessions"],
        prefix=lambda wc: f"{OUTPUT_ROOT}/{wc.db}/supplemental",
        taxdump_dir=TAXDUMP_DIR,
    threads: lambda wc: res(wc.db, "fetch_taxonomy", "threads", 2)
    resources:
        mem_gb=lambda wc: res(wc.db, "fetch_taxonomy", "mem_gb", 4),
        time=lambda wc: res(wc.db, "fetch_taxonomy", "time", 120),
        partition=lambda wc: res(wc.db, "fetch_taxonomy", "partition", "compute"),
    conda: "workflow/env/ncbi_taxonkit.yaml"
    shell:
        """
        mkdir -p $(dirname {output.tax})
        python workflow/scripts/ncbi_to_gtdb_taxonomy.py \
            --accessions {params.accessions} \
            --out-prefix {params.prefix} \
            --taxdump-dir {params.taxdump_dir}
        """


# ---- GTDB downloads (only when gtdb_release is configured) ------------------

rule get_gtdb_bac_taxonomy:
    output: f"{OUTPUT_ROOT}/{{db}}/gtdb/bac120_taxonomy.tsv"
    params: url=lambda wc: gtdb_url(wc.db, "bac_tax_url")
    threads: 2
    resources: mem_gb=4, time=120, partition="compute"
    shell:
        """
        mkdir -p $(dirname {output})
        curl -fsSL -o {output}.gz "{params.url}"
        gunzip -f {output}.gz
        """


rule get_gtdb_arc_taxonomy:
    output: f"{OUTPUT_ROOT}/{{db}}/gtdb/ar53_taxonomy.tsv"
    params: url=lambda wc: gtdb_url(wc.db, "arc_tax_url")
    threads: 2
    resources: mem_gb=4, time=120, partition="compute"
    shell:
        """
        mkdir -p $(dirname {output})
        curl -fsSL -o {output}.gz "{params.url}"
        gunzip -f {output}.gz
        """


rule get_gtdb_bac_metadata:
    output: f"{OUTPUT_ROOT}/{{db}}/gtdb/bac120_metadata.tsv"
    params: url=lambda wc: gtdb_url(wc.db, "bac_meta_url")
    threads: 2
    resources: mem_gb=4, time=120, partition="compute"
    shell:
        """
        mkdir -p $(dirname {output})
        curl -fsSL -o {output}.gz "{params.url}"
        gunzip -f {output}.gz
        """


rule get_gtdb_arc_metadata:
    output: f"{OUTPUT_ROOT}/{{db}}/gtdb/ar53_metadata.tsv"
    params: url=lambda wc: gtdb_url(wc.db, "arc_meta_url")
    threads: 2
    resources: mem_gb=4, time=120, partition="compute"
    shell:
        """
        mkdir -p $(dirname {output})
        curl -fsSL -o {output}.gz "{params.url}"
        gunzip -f {output}.gz
        """


rule make_merged_metadata:
    output: f"{OUTPUT_ROOT}/{{db}}/gtdb/merged_metadata.tsv"
    input:
        bac=f"{OUTPUT_ROOT}/{{db}}/gtdb/bac120_metadata.tsv",
        arc=f"{OUTPUT_ROOT}/{{db}}/gtdb/ar53_metadata.tsv",
    threads: 2
    resources: mem_gb=4, time=60, partition="compute"
    shell:
        """
        cat {input.bac} > {output}
        grep -v '^accession' {input.arc} >> {output}
        """


rule download_gtdb_genomes_tar:
    output: f"{OUTPUT_ROOT}/{{db}}/gtdb/gtdb_genomes_reps.tar.gz"
    params: url=lambda wc: gtdb_url(wc.db, "genomes_url")
    threads: 4
    resources: mem_gb=8, time=720, partition="compute"
    shell:
        """
        mkdir -p $(dirname {output})
        curl -fsSL --retry 5 --retry-delay 30 -o {output} "{params.url}"
        """


# ---- Merged taxonomy + per-accession downloads ------------------------------

rule make_merged_taxonomy:
    output: f"{OUTPUT_ROOT}/{{db}}/merged_taxonomy.tsv"
    input: unpack(merged_tax_inputs)
    threads: 2
    resources: mem_gb=4, time=60, partition="compute"
    run:
        with open(output[0], "w") as out_fh:
            if "bac" in input.keys():
                with open(input.bac) as fh:
                    for line in fh:
                        out_fh.write(line)
            if "arc" in input.keys():
                with open(input.arc) as fh:
                    next(fh, None)  # skip header if present
                    for line in fh:
                        out_fh.write(line)
            with open(input.supp) as fh:
                for line in fh:
                    out_fh.write(line)


rule download_supplemental_genome:
    output: f"{OUTPUT_ROOT}/{{db}}/downloads/{{accession}}.fna.gz"
    input:
        urls=f"{OUTPUT_ROOT}/{{db}}/supplemental.urls.tsv",
    threads: 1
    resources: mem_gb=2, time=120, partition="compute"
    shell:
        """
        mkdir -p $(dirname {output})
        url=$(awk -F'\\t' -v acc='{wildcards.accession}' '$1==acc {{print $2}}' {input.urls})
        if [ -z "$url" ]; then
            echo "Error: no URL for accession {wildcards.accession} in {input.urls}" >&2
            exit 1
        fi
        curl -fsSL --retry 5 --retry-delay 15 -o {output} "$url"
        """


# ---- Collect, convert, build ------------------------------------------------

def collect_inputs(wildcards):
    db = wildcards.db
    accs = accessions_of(db)
    inputs = {
        "supplemental": expand(
            f"{OUTPUT_ROOT}/{{db}}/downloads/{{acc}}.fna.gz",
            db=[db], acc=accs,
        ),
    }
    if has_gtdb(db):
        inputs["gtdb_tar"] = f"{OUTPUT_ROOT}/{db}/gtdb/gtdb_genomes_reps.tar.gz"
    return inputs


rule collect_genomes:
    output:
        flag=f"{OUTPUT_ROOT}/{{db}}/input_genomes/.collected",
        dir=directory(f"{OUTPUT_ROOT}/{{db}}/input_genomes"),
    input: unpack(collect_inputs)
    params:
        gtdb_tar=lambda wc: (f"{OUTPUT_ROOT}/{wc.db}/gtdb/gtdb_genomes_reps.tar.gz"
                             if has_gtdb(wc.db) else ""),
    threads: 4
    resources: mem_gb=8, time=240, partition="compute"
    shell:
        """
        mkdir -p {output.dir}
        for f in {input.supplemental}; do
            cp -n "$f" {output.dir}/
        done
        if [ -n "{params.gtdb_tar}" ] && [ -s "{params.gtdb_tar}" ]; then
            tar -xzf "{params.gtdb_tar}" -C {output.dir}
            find {output.dir} -mindepth 2 -type f -name '*.fna.gz' -exec mv -n {{}} {output.dir}/ \\;
            find {output.dir} -mindepth 1 -type d -empty -delete
        fi
        touch {output.flag}
        """


rule convert_taxonomy:
    output:
        nodes=f"{OUTPUT_ROOT}/{{db}}/taxonomy/nodes.dmp",
        names=f"{OUTPUT_ROOT}/{{db}}/taxonomy/names.dmp",
        kraken_dir=directory(f"{OUTPUT_ROOT}/{{db}}/kraken_input"),
        conv=f"{OUTPUT_ROOT}/{{db}}/taxonomy/conversion.tsv",
    input:
        merged=f"{OUTPUT_ROOT}/{{db}}/merged_taxonomy.tsv",
        flag=f"{OUTPUT_ROOT}/{{db}}/input_genomes/.collected",
        genomes_dir=f"{OUTPUT_ROOT}/{{db}}/input_genomes",
    threads: 4
    resources: mem_gb=16, time=240, partition="compute"
    conda: "workflow/env/kraken2.yaml"
    shell:
        """
        mkdir -p $(dirname {output.nodes})
        rm -rf {output.kraken_dir}
        python workflow/scripts/tax_from_gtdb.py \
            --gtdb {input.merged} \
            --assemblies {input.genomes_dir} \
            --nodes {output.nodes} \
            --names {output.names} \
            --conversion {output.conv} \
            --kraken_dir {output.kraken_dir}
        """


rule populate_library:
    output:
        flag=f"{OUTPUT_ROOT}/{{db}}/.library_populated",
    input:
        kraken_dir=f"{OUTPUT_ROOT}/{{db}}/kraken_input",
        nodes=f"{OUTPUT_ROOT}/{{db}}/taxonomy/nodes.dmp",
        names=f"{OUTPUT_ROOT}/{{db}}/taxonomy/names.dmp",
    params:
        db_dir=lambda wc: f"{OUTPUT_ROOT}/{wc.db}",
    threads: lambda wc: res(wc.db, "add_to_library", "threads", 16)
    resources:
        mem_gb=lambda wc: res(wc.db, "add_to_library", "mem_gb", 24),
        time=lambda wc: res(wc.db, "add_to_library", "time", 1440),
        partition=lambda wc: res(wc.db, "add_to_library", "partition", "compute"),
    conda: "workflow/env/kraken2.yaml"
    shell:
        """
        # nodes.dmp/names.dmp are already at {params.db_dir}/taxonomy/ from convert_taxonomy
        rm -rf {params.db_dir}/library
        shopt -s nullglob
        for fa in {input.kraken_dir}/*.fa {input.kraken_dir}/*.fna; do
            kraken2-build --add-to-library "$fa" --db {params.db_dir} --threads {threads} --no-masking
        done
        touch {output.flag}
        """


rule build_kraken2:
    output:
        hash=f"{OUTPUT_ROOT}/{{db}}/hash.k2d",
        opts=f"{OUTPUT_ROOT}/{{db}}/opts.k2d",
        taxo=f"{OUTPUT_ROOT}/{{db}}/taxo.k2d",
    input:
        flag=f"{OUTPUT_ROOT}/{{db}}/.library_populated",
    params:
        db_dir=lambda wc: f"{OUTPUT_ROOT}/{wc.db}",
        kmer=lambda wc: spec_of(wc.db).get("build_params", {}).get("kmer_len", 35),
        mink=lambda wc: spec_of(wc.db).get("build_params", {}).get("minimizer_len", 31),
        sp=lambda wc: spec_of(wc.db).get("build_params", {}).get("minimizer_spaces", 7),
    threads: lambda wc: res(wc.db, "build", "threads", 64)
    resources:
        mem_gb=lambda wc: res(wc.db, "build", "mem_gb", 1600),
        time=lambda wc: res(wc.db, "build", "time", 4320),
        partition=lambda wc: res(wc.db, "build", "partition", "hugemem"),
    conda: "workflow/env/kraken2.yaml"
    shell:
        """
        kraken2-build --build \
            --db {params.db_dir} \
            --threads {threads} \
            --kmer-len {params.kmer} \
            --minimizer-len {params.mink} \
            --minimizer-spaces {params.sp}
        """


# ---- Manifest --------------------------------------------------------------

def manifest_inputs(wildcards):
    db = wildcards.db
    base = f"{OUTPUT_ROOT}/{db}"
    inputs = {
        "hash": f"{base}/hash.k2d",
        "supp_tax": f"{base}/supplemental.taxonomy.tsv",
        "supp_urls": f"{base}/supplemental.urls.tsv",
        "merged_tax": f"{base}/merged_taxonomy.tsv",
    }
    if has_gtdb(db):
        inputs["bac_tax"] = f"{base}/gtdb/bac120_taxonomy.tsv"
        inputs["arc_tax"] = f"{base}/gtdb/ar53_taxonomy.tsv"
        inputs["bac_meta"] = f"{base}/gtdb/bac120_metadata.tsv"
        inputs["arc_meta"] = f"{base}/gtdb/ar53_metadata.tsv"
        inputs["gtdb_tar"] = f"{base}/gtdb/gtdb_genomes_reps.tar.gz"
    return inputs


rule write_manifest:
    output: f"{OUTPUT_ROOT}/{{db}}/MANIFEST.{{db}}.json"
    input: unpack(manifest_inputs)
    params:
        db_dir=lambda wc: f"{OUTPUT_ROOT}/{wc.db}",
        name=lambda wc: wc.db,
        label=lambda wc: spec_of(wc.db).get("label", DB_BY_DIR[wc.db]),
        release=lambda wc: spec_of(wc.db).get("gtdb_release", "") or "",
        accessions_file=lambda wc: spec_of(wc.db)["supplemental_accessions"],
        kmer=lambda wc: spec_of(wc.db).get("build_params", {}).get("kmer_len", 35),
        mink=lambda wc: spec_of(wc.db).get("build_params", {}).get("minimizer_len", 31),
        sp=lambda wc: spec_of(wc.db).get("build_params", {}).get("minimizer_spaces", 7),
        opt_gtdb_args=lambda wc: (
            f"--gtdb-bac-tax {OUTPUT_ROOT}/{wc.db}/gtdb/bac120_taxonomy.tsv "
            f"--gtdb-arc-tax {OUTPUT_ROOT}/{wc.db}/gtdb/ar53_taxonomy.tsv "
            f"--gtdb-bac-meta {OUTPUT_ROOT}/{wc.db}/gtdb/bac120_metadata.tsv "
            f"--gtdb-arc-meta {OUTPUT_ROOT}/{wc.db}/gtdb/ar53_metadata.tsv "
            f"--gtdb-genomes-tar {OUTPUT_ROOT}/{wc.db}/gtdb/gtdb_genomes_reps.tar.gz"
        ) if has_gtdb(wc.db) else "",
    threads: 1
    resources: mem_gb=2, time=15, partition="compute"
    conda: "workflow/env/kraken2.yaml"
    shell:
        """
        python workflow/scripts/write_manifest.py \
            --db-dir {params.db_dir} \
            --name {params.name} \
            --label {params.label} \
            --gtdb-release "{params.release}" \
            --accessions-file {params.accessions_file} \
            --supplemental-taxonomy {input.supp_tax} \
            --supplemental-urls {input.supp_urls} \
            --merged-taxonomy {input.merged_tax} \
            --kmer-len {params.kmer} \
            --minimizer-len {params.mink} \
            --minimizer-spaces {params.sp} \
            {params.opt_gtdb_args}
        """
