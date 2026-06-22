#!/usr/bin/env python3
"""Convert a list of NCBI accessions into GTDB-style taxonomy + a URL manifest.

For each accession (GCF_* or GCA_*) it:
  1. Calls `datasets summary genome accession ...` to get tax_id and the genomic FNA URL
  2. Resolves the NCBI taxid to a GTDB-style lineage via `taxonkit lineage | taxonkit reformat -P`
  3. Emits two TSVs:
        <out_prefix>.taxonomy.tsv   GTDB-style taxonomy (matches tax_from_gtdb.py input)
        <out_prefix>.urls.tsv       accession<TAB>https://ftp.ncbi.nlm.nih.gov/.../genomic.fna.gz

Prefix rules: GCF_* -> RS_, GCA_* -> GB_ (matches tax_from_gtdb.py's expectation).

Requires: ncbi-datasets-cli >=16, taxonkit >=0.13. Optionally NCBI_API_KEY in env.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path


_ACC_RE = re.compile(r"^(GC[AF])_\d{9}\.\d+$")
_BATCH = 50


def read_accessions(path: Path) -> list[str]:
    accs: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        token = line.split()[0]
        if not _ACC_RE.match(token):
            sys.exit(f"Error: {path}: invalid accession '{token}' (need GCF_NNNNNNNNN.V or GCA_NNNNNNNNN.V)")
        accs.append(token)
    seen = set()
    out = []
    for a in accs:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def run_with_retry(cmd: list[str], stdin: str | None = None, attempts: int = 4) -> str:
    last_err = None
    for i in range(attempts):
        try:
            proc = subprocess.run(
                cmd,
                input=stdin,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            return proc.stdout
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            last_err = exc
            wait = 2 ** i
            print(f"  retry {i + 1}/{attempts} after {wait}s ({type(exc).__name__})", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"command failed after {attempts} attempts: {cmd!r}\nlast error: {last_err}")


def fetch_summaries(accessions: list[str]) -> dict[str, dict]:
    """Batch-call `datasets summary genome accession` and return a dict keyed by accession."""
    summaries: dict[str, dict] = {}
    for i in range(0, len(accessions), _BATCH):
        batch = accessions[i:i + _BATCH]
        print(f"datasets summary batch {i // _BATCH + 1}: {len(batch)} accessions", file=sys.stderr)
        cmd = ["datasets", "summary", "genome", "accession", *batch, "--as-json-lines"]
        out = run_with_retry(cmd)
        for line in out.splitlines():
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            acc = rec.get("accession")
            if acc:
                summaries[acc] = rec
    missing = [a for a in accessions if a not in summaries]
    if missing:
        sys.exit(f"Error: datasets returned no summary for: {missing}")
    return summaries


def construct_ftp_url(accession: str, asm_name: str) -> str:
    """Build the canonical NCBI FTP URL for a genome assembly's .fna.gz file."""
    prefix, digits_version = accession.split("_", 1)  # ('GCF', '016772045.2')
    digits = digits_version.split(".", 1)[0]
    parts = [digits[0:3], digits[3:6], digits[6:9]]
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{parts[0]}/{parts[1]}/{parts[2]}"
    folder = f"{accession}_{asm_name}"
    fname = f"{accession}_{asm_name}_genomic.fna.gz"
    return f"{base}/{folder}/{fname}"


def extract_taxid_and_url(rec: dict, accession: str) -> tuple[str, str]:
    taxid = None
    org = rec.get("organism") or {}
    taxid = org.get("tax_id") or rec.get("tax_id")
    if not taxid:
        sys.exit(f"Error: no tax_id for {accession} in datasets summary")

    asm_info = rec.get("assembly_info") or {}
    asm_name = asm_info.get("assembly_name") or rec.get("assembly_name")
    if not asm_name:
        sys.exit(f"Error: no assembly_name for {accession} in datasets summary")
    asm_name = asm_name.replace(" ", "_")

    url = construct_ftp_url(accession, asm_name)
    return str(taxid), url


def taxonkit_lineages(taxids: list[str], data_dir: str | None = None) -> dict[str, str]:
    """Map taxid -> GTDB-style 7-rank lineage string."""
    stdin = "\n".join(taxids) + "\n"
    fmt = "d__{k};p__{p};c__{c};o__{o};f__{f};g__{g};s__{s}"
    extra = ["--data-dir", data_dir] if data_dir else []
    cmd_lineage = ["taxonkit", *extra, "lineage"]
    cmd_reformat = ["taxonkit", *extra, "reformat", "-P", "-f", fmt, "--miss-rank-repl", "unclassified", "-I", "1"]

    p1 = subprocess.Popen(cmd_lineage, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    p2 = subprocess.Popen(cmd_reformat, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    p1.stdout.close()
    p1.stdin.write(stdin)
    p1.stdin.close()
    out, err = p2.communicate()
    rc1 = p1.wait()
    rc2 = p2.returncode
    if rc1 != 0 or rc2 != 0:
        sys.exit(f"taxonkit failed (lineage rc={rc1}, reformat rc={rc2}):\n{err}")

    result: dict[str, str] = {}
    for line in out.splitlines():
        cols = line.split("\t")
        if len(cols) < 3:
            continue
        taxid = cols[0].strip()
        lineage = cols[-1].strip()
        if not lineage.startswith("d__"):
            sys.exit(f"taxonkit reformat returned non-conforming lineage for taxid {taxid}: '{lineage}'")
        result[taxid] = lineage
    missing = [t for t in taxids if t not in result]
    if missing:
        sys.exit(f"taxonkit returned no lineage for taxids: {missing}")
    return result


def prefix_for(accession: str) -> str:
    return "RS_" if accession.startswith("GCF_") else "GB_"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", type=Path, required=True)
    ap.add_argument("--out-prefix", type=Path, required=True,
                    help="Output prefix; writes <prefix>.taxonomy.tsv and <prefix>.urls.tsv")
    ap.add_argument("--taxdump-dir", type=str, default=None,
                    help="Directory containing nodes.dmp/names.dmp/merged.dmp/delnodes.dmp "
                         "(passed as --data-dir to taxonkit). Defaults to ~/.taxonkit/")
    args = ap.parse_args()

    if not os.environ.get("NCBI_API_KEY"):
        print("Note: NCBI_API_KEY not set; rate limit is ~3 req/s anonymous.", file=sys.stderr)

    accs = read_accessions(args.accessions)
    print(f"Processing {len(accs)} accessions", file=sys.stderr)

    summaries = fetch_summaries(accs)

    taxid_per_acc: dict[str, str] = {}
    url_per_acc: dict[str, str] = {}
    for acc in accs:
        taxid, url = extract_taxid_and_url(summaries[acc], acc)
        taxid_per_acc[acc] = taxid
        url_per_acc[acc] = url

    unique_taxids = sorted(set(taxid_per_acc.values()))
    lineages = taxonkit_lineages(unique_taxids, data_dir=args.taxdump_dir)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tax_path = args.out_prefix.with_suffix(args.out_prefix.suffix + ".taxonomy.tsv") \
        if args.out_prefix.suffix else Path(str(args.out_prefix) + ".taxonomy.tsv")
    urls_path = args.out_prefix.with_suffix(args.out_prefix.suffix + ".urls.tsv") \
        if args.out_prefix.suffix else Path(str(args.out_prefix) + ".urls.tsv")

    with tax_path.open("w") as fh:
        for acc in accs:
            fh.write(f"{prefix_for(acc)}{acc}\t{lineages[taxid_per_acc[acc]]}\n")

    with urls_path.open("w") as fh:
        for acc in accs:
            fh.write(f"{acc}\t{url_per_acc[acc]}\n")

    print(f"wrote {tax_path}", file=sys.stderr)
    print(f"wrote {urls_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
