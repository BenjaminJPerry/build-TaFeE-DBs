# build-GTDB-DBs

Reproducible Snakemake workflow that builds **Kraken2 indices** from any user-supplied list of NCBI reference genome accessions, optionally combined with a specified GTDB release. The workflow auto-derives all download URLs and taxonomy, names the output directory from the input provenance, and writes a full provenance manifest beside the index.

The same Snakefile produces (a) hosts-only databases for read decontamination, (b) GTDB-based prokaryote databases with optional supplemental eukaryotes, or (c) any other combination — the only inputs are an accessions text file and (optionally) a GTDB release version.

## Quickstart

```bash
module load snakemake
cd /mnt/gpfs/persist/projects/2022-bjp-gtdb/build-GTDB-DBs
snakemake --configfile config/test_config.yaml --profile eRI
```

A typical hosts-only run takes ~3 hours wall, peaks at ~20 GB RAM, and emits a ~12 GB `hash.k2d`. SLURM jobs are submitted by the controller on the login node.

## Database naming

Each database produces a directory whose name encodes its provenance:

```
kraken2-gtdb-r{release}-{label}-{shortsha}/      # GTDB-based builds
kraken2-{label}-{shortsha}/                       # no-GTDB (e.g. hosts-only)
```

- `release`  — value of `gtdb_release` in the config (e.g. `220.0`)
- `label`    — free-form slug from the per-DB config (e.g. `hosts`, `tickGBS`)
- `shortsha` — first 8 hex chars of `sha256` over the sorted, deduped accession list

Different accession lists against the same release produce distinct, non-colliding directories. A sidecar `MANIFEST.<db>.json` (e.g. `MANIFEST.kraken2-hosts-87fed59e.json`) next to `hash.k2d` records the full provenance (URLs, file sha256s, kraken2 / datasets / taxonkit versions, git commit, build params, accession list). The filename embeds the database name so the manifest is self-identifying when copied out of its directory.

The optional `name_override:` per-DB key bypasses auto-naming when an exact legacy path must be preserved.

## Rule graph

The DAG branches at DAG-construction time on whether a database has `gtdb_release` set. GTDB-only rules are silently skipped when it is absent.

```
                ┌─────────────────────┐
                │  download_taxdump   │  shared NCBI taxdump.tar.gz
                └──────────┬──────────┘
                           │ nodes/names/merged/delnodes.dmp
                           ▼
   ┌───────────────────────────────────┐  per-accession metadata
   │  fetch_supplemental_taxonomy      │  via `datasets summary` + `taxonkit`
   └───────────────────┬───────────────┘
                       │ supplemental.taxonomy.tsv + supplemental.urls.tsv
                       │
   (optional, when gtdb_release is set)
   ┌──────────────────┴───────────────┐
   │ get_gtdb_bac_taxonomy             │
   │ get_gtdb_arc_taxonomy             │
   │ get_gtdb_bac_metadata             │
   │ get_gtdb_arc_metadata             │
   │ download_gtdb_genomes_tar         │
   │ make_merged_metadata              │
   └──────────────────┬───────────────┘
                       │
                       ▼
              ┌─────────────────────┐         ┌───────────────────────────┐
              │ make_merged_taxonomy│         │ download_supplemental_    │
              └──────────┬──────────┘         │    genome (1 per accession)│
                         │                    └────────────┬──────────────┘
                         │                                 │
                         └───────────┐  ┌──────────────────┘
                                     ▼  ▼
                              ┌────────────────┐
                              │ collect_genomes│  populates input_genomes/
                              └────────┬───────┘
                                       ▼
                             ┌────────────────────┐
                             │ convert_taxonomy   │  tax_from_gtdb.py:
                             │                    │  -> nodes.dmp, names.dmp,
                             │                    │     kraken_input/*.fa
                             └────────┬───────────┘
                                      ▼
                             ┌────────────────────┐
                             │ populate_library   │  kraken2-build --add-to-library
                             └────────┬───────────┘
                                      ▼
                             ┌────────────────────┐
                             │ build_kraken2      │  kraken2-build --build
                             └────────┬───────────┘
                                      │  hash.k2d, opts.k2d, taxo.k2d
                                      ▼
                             ┌────────────────────┐
                             │ write_manifest     │  MANIFEST.<db>.json
                             └────────────────────┘
```

### Per-rule reference

Paths are relative to `{output_root}/{db}/` (where `{db}` is the auto-named directory) unless noted.

#### `download_taxdump`
Downloads NCBI's taxdump (`nodes.dmp`, `names.dmp`, `merged.dmp`, `delnodes.dmp`) to a shared location used by every database in the same `output_root`. Required by `taxonkit` for taxid → lineage resolution.

- **Inputs**: none (fetches from `ftp.ncbi.nlm.nih.gov`)
- **Outputs**: `{output_root}/_taxdump/{nodes,names,merged,delnodes}.dmp`
- **Env**: shell only (`curl`, `tar`)

#### `fetch_supplemental_taxonomy`
For every accession in the per-DB list, queries `datasets summary genome accession` to get NCBI taxid and assembly name, then converts each taxid into a GTDB-style 7-rank lineage via `taxonkit reformat -P`. Emits one TSV in the format `tax_from_gtdb.py` consumes, plus a parallel TSV of FTP URLs for the genome downloads.

- **Inputs**: `{output_root}/_taxdump/nodes.dmp`, accessions file from config
- **Outputs**: `supplemental.taxonomy.tsv`, `supplemental.urls.tsv`
- **Env**: `workflow/env/ncbi_taxonkit.yaml` (`ncbi-datasets-cli`, `taxonkit`)
- **Prefix rule**: `GCF_*` → `RS_`, `GCA_*` → `GB_` (matches `tax_from_gtdb.py`)

#### `get_gtdb_bac_taxonomy` / `get_gtdb_arc_taxonomy` *(GTDB-only)*
Downloads the per-release bacterial and archaeal taxonomy TSVs from `data.ace.uq.edu.au`. URLs derived from `gtdb_release` by `workflow/scripts/resolve_gtdb_urls.py`.

- **Outputs**: `gtdb/bac120_taxonomy.tsv`, `gtdb/ar53_taxonomy.tsv`
- **Env**: shell only

#### `get_gtdb_bac_metadata` / `get_gtdb_arc_metadata` *(GTDB-only)*
Same as above for the metadata TSVs (used by `make_merged_metadata` and downstream Struo2-style tooling).

- **Outputs**: `gtdb/bac120_metadata.tsv`, `gtdb/ar53_metadata.tsv`

#### `make_merged_metadata` *(GTDB-only)*
Concatenates bac + arc metadata into a single file. Provided for downstream tools that expect a unified metadata table; not required by the kraken2 build itself.

- **Inputs**: `gtdb/bac120_metadata.tsv`, `gtdb/ar53_metadata.tsv`
- **Outputs**: `gtdb/merged_metadata.tsv`

#### `make_merged_taxonomy`
Concatenates GTDB bac + arc taxonomy (when `gtdb_release` is set) with the supplemental taxonomy from `fetch_supplemental_taxonomy`. For hosts-only builds this is just the supplemental file.

- **Inputs**: `supplemental.taxonomy.tsv` (always), `gtdb/bac120_taxonomy.tsv` + `gtdb/ar53_taxonomy.tsv` (conditional)
- **Outputs**: `merged_taxonomy.tsv`

#### `download_gtdb_genomes_tar` *(GTDB-only)*
Downloads the GTDB representative genomes tarball for the configured release. This is the large file (~80 GB for r220).

- **Outputs**: `gtdb/gtdb_genomes_reps.tar.gz`

#### `download_supplemental_genome`
One job per accession. Reads the FTP URL from `supplemental.urls.tsv` and downloads the `_genomic.fna.gz` file.

- **Inputs**: `supplemental.urls.tsv`
- **Outputs**: `downloads/{accession}.fna.gz` (one per accession)

#### `collect_genomes`
Gathers all genome `.fna.gz` files into a single flat directory `input_genomes/`. For GTDB-enabled builds, extracts the GTDB reps tarball and flattens nested directories.

- **Inputs**: all `downloads/*.fna.gz` (via `expand` over the accession list) + the GTDB tar (conditional)
- **Outputs**: `input_genomes/` (directory) + `input_genomes/.collected` flag

#### `convert_taxonomy`
Runs the existing `workflow/scripts/tax_from_gtdb.py` to turn the merged GTDB-style taxonomy + the flat `input_genomes/` directory into:
- NCBI-style `nodes.dmp` and `names.dmp`
- a `kraken_input/` directory of per-accession FASTAs whose headers are rewritten to include `|kraken:taxid|<id>` markers
- a `conversion.tsv` mapping contigs to taxids

- **Inputs**: `merged_taxonomy.tsv`, `input_genomes/`, `input_genomes/.collected`
- **Outputs**: `taxonomy/nodes.dmp`, `taxonomy/names.dmp`, `taxonomy/conversion.tsv`, `kraken_input/`
- **Env**: `workflow/env/kraken2.yaml` (needs `biopython`)

#### `populate_library`
Iterates over `kraken_input/*.fa` and calls `kraken2-build --add-to-library` for each. The taxonomy files are already in `taxonomy/` from the previous rule, so kraken2-build finds them automatically.

- **Inputs**: `kraken_input/`, `taxonomy/nodes.dmp`, `taxonomy/names.dmp`
- **Outputs**: `.library_populated` flag (the cumulative state is in `library/`)
- **Env**: `workflow/env/kraken2.yaml`

#### `build_kraken2`
The expensive step. Runs `kraken2-build --build` with the configured `kmer_len`, `minimizer_len`, and `minimizer_spaces`. Produces the final hash table.

- **Inputs**: `.library_populated`
- **Outputs**: `hash.k2d`, `opts.k2d`, `taxo.k2d`
- **Env**: `workflow/env/kraken2.yaml`
- **Resources**: from per-DB `resources.build` in the config

#### `write_manifest`
Final rule. Writes `MANIFEST.<db>.json` next to `hash.k2d` with the full provenance: gtdb release (if any), URLs + sha256 of every input file, the verbatim accession list, kraken2 / taxonkit / datasets versions, the build host, git commit, and kraken2 build params. The filename embeds the auto-computed db name (e.g. `MANIFEST.kraken2-hosts-87fed59e.json`) so the manifest is self-identifying.

- **Inputs**: `hash.k2d`, all upstream TSVs and (for GTDB builds) the GTDB sources
- **Outputs**: `MANIFEST.<db>.json`
- **Env**: `workflow/env/kraken2.yaml` (provides `kraken2`, `ncbi-datasets-cli`, and `taxonkit` so the manifest records all four tool versions)

#### `targets`
The default rule. Requests `hash.k2d` + `MANIFEST.<db>.json` for every database defined in the config.

## Adding a new database

The workflow is designed so that adding a database requires **no Snakefile changes** — only a new accession list and a config block.

### 1. Create the accession list

Plain text, one accession per line. `#` introduces a comment; blank lines are ignored. Optional second column is a free-form label (preserved in the manifest, ignored by the build). Accessions **must** include the version suffix.

```text
# resources/genomes/my_taxa.txt
# Mixed RefSeq and GenBank accessions for the FOO project.
GCF_016772045.2  sheep          # mammal host
GCF_002263795.3  cow            # mammal host
GCA_002087855.3  Entodinium     # ciliate
```

Accepted formats: `GCF_NNNNNNNNN.V` (RefSeq) and `GCA_NNNNNNNNN.V` (GenBank).

### 2. Add a `databases:` block

Edit `config/config.yaml` (or create a new config file). Each entry under `databases:` is one DB:

```yaml
output_root: GTDB                       # parent directory for all DBs

databases:

  my_db:                                # arbitrary key, used only internally
    label: foo-project                  # appears in the auto-generated directory name
    gtdb_release: "220.0"               # omit for a no-GTDB build
    supplemental_accessions: resources/genomes/my_taxa.txt
    # name_override: legacy-name        # optional: skip auto-naming
    build_params:
      kmer_len: 35
      minimizer_len: 31
      minimizer_spaces: 7
    resources:
      fetch_taxonomy: {threads: 2,  mem_gb: 4,  time: 60,   partition: compute}
      add_to_library: {threads: 8,  mem_gb: 16, time: 120,  partition: compute}
      build:          {threads: 32, mem_gb: 40, time: 1440, partition: "compute,hugemem"}
```

The resulting directory will be:
- `GTDB/kraken2-gtdb-r220.0-foo-project-<8-hex-shortsha>/` (with GTDB)
- `GTDB/kraken2-foo-project-<8-hex-shortsha>/` (no `gtdb_release`)

The `partition` value can be a comma-separated list (e.g. `"compute,hugemem"`) — SLURM will land the job on whichever has resources first.

### 3. Dry-run to validate

```bash
module load snakemake
snakemake --configfile config/config.yaml --profile eRI -n
```

Inspect the printed `db '<key>' -> <output_root>/<computed-name>` line, the job count, and the resources for `build_kraken2`. Any errors here (missing accessions file, invalid accession format) surface before any jobs are submitted.

### 4. Pre-build conda envs *(first run only)*

The workflow uses two conda envs (`workflow/env/ncbi_taxonkit.yaml`, `workflow/env/kraken2.yaml`). Pre-build them once to avoid every SLURM job retrying env creation:

```bash
snakemake --configfile config/config.yaml --use-conda --conda-create-envs-only --cores 4
```

This takes ~10 minutes total. The built envs land in `.snakemake/conda/` and are reused across runs.

### 5. Launch

```bash
nohup snakemake --configfile config/config.yaml --profile eRI \
    > logs/snakemake-$(date +%Y%m%d-%H%M%S).log 2>&1 &
echo "controller PID: $!"
```

The Snakemake controller stays on the login node and submits each rule as a SLURM job via the `eRI/` profile. Job state is visible via `squeue -u $USER`; per-rule logs land in `logs/<rule>/<rule>-<wildcards>-<jobid>.out`.

### 6. Watch progress

```bash
# Snakemake controller log (overall DAG progress)
tail -f $(/bin/ls -t logs/snakemake-*.log | head -1)

# SLURM queue
squeue -u $USER -o "%i %P %T %r %j %M %l"

# Per-rule outputs
ls -lh GTDB/kraken2-*/hash.k2d 2>/dev/null
```

When the controller finishes you'll see `N of N steps (100%) done`. The index is then ready at `{output_root}/{computed-name}/hash.k2d`, with full provenance in `MANIFEST.{computed-name}.json` alongside.

## Testing on a small DB (hosts only)

`config/test_config.yaml` defines a hosts-only build into `test_run/` (so it can't touch any existing `GTDB/` artifacts):

```bash
module load snakemake
snakemake --configfile config/test_config.yaml --profile eRI -n     # dry-run
nohup snakemake --configfile config/test_config.yaml --profile eRI \
    > logs/snakemake-hosts-$(date +%Y%m%d-%H%M%S).log 2>&1 &
```

Benchmark (`benchmark/build_kraken2_hosts.txt`): 5 mammalian hosts, peak RSS 19 GB, wall 2h47m, 12 GB `hash.k2d`.

## Layout

```
build-GTDB-DBs/
├── Snakefile
├── config/
│   ├── config.yaml             # production DBs (e.g. supp220, hosts)
│   └── test_config.yaml        # hosts-only build into test_run/
├── eRI/
│   ├── config.yaml             # SLURM profile (snakemake 9, cluster-generic plugin)
│   └── ...
├── resources/
│   ├── genomes/
│   │   ├── hosts_accessions.txt
│   │   └── GTDB-220.0_supplemental_accessions.txt
│   └── host_taxa.tsv           # legacy reference, no longer consumed
├── workflow/
│   ├── env/
│   │   ├── ncbi_taxonkit.yaml
│   │   └── kraken2.yaml
│   └── scripts/
│       ├── compute_db_name.py          # accessions -> shortsha + directory name
│       ├── resolve_gtdb_urls.py        # release version -> 5 canonical URLs
│       ├── ncbi_to_gtdb_taxonomy.py    # datasets + taxonkit -> taxonomy + urls TSVs
│       ├── tax_from_gtdb.py            # GTDB-style TSV -> nodes/names/kraken_dir
│       ├── write_manifest.py           # final MANIFEST.<db>.json emitter
│       └── status.py                   # SLURM job-status helper for the profile
├── GTDB/                       # default output_root for production DBs
└── test_run/                   # default output_root for the test config
```

## Cluster profile notes

`eRI/config.yaml` uses the Snakemake 9 `cluster-generic` executor plugin. Resources flow through as:
- `--cpus-per-task={threads}`
- `--mem={resources.mem_gb}G`
- `--time={resources.time}` (minutes)
- `--partition={resources.partition}`
- `--account={resources.account}` (defaults to `2023-mbie-rumen-gbs`)

The legacy Snakemake 7 profile (using the `cluster:` directive) is preserved at `eRI-v7-backup.yaml` (kept outside `eRI/` because Snakemake reads every file in the profile directory).

## Future work

Out of scope for this iteration:
- humann3 functional databases (still manual)
- SILVA bowtie indices for Kneaddata (still manual)
- Prebuilt NCBI nt download (kept as a separate dedicated step)
