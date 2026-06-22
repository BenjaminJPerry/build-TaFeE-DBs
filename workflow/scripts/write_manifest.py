#!/usr/bin/env python3
"""Write MANIFEST.<name>.json next to hash.k2d capturing full build provenance."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def get_version(cmd: list[str]) -> str:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20, check=False).stdout.strip()
        return out.splitlines()[0] if out else "unknown"
    except Exception as exc:
        return f"unknown ({exc})"


def git_commit(repo: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo, capture_output=True, text=True, timeout=10, check=False,
        )
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def file_record(path: Path, key: str) -> dict | None:
    if path is None or not path.exists():
        return None
    return {"path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-dir", type=Path, required=True, help="Directory containing hash.k2d")
    ap.add_argument("--name", type=str, required=True)
    ap.add_argument("--label", type=str, required=True)
    ap.add_argument("--gtdb-release", type=str, default="")
    ap.add_argument("--accessions-file", type=Path, required=True)
    ap.add_argument("--supplemental-taxonomy", type=Path, required=True)
    ap.add_argument("--supplemental-urls", type=Path, required=True)
    ap.add_argument("--merged-taxonomy", type=Path, required=True)
    ap.add_argument("--gtdb-bac-tax", type=Path, default=None)
    ap.add_argument("--gtdb-arc-tax", type=Path, default=None)
    ap.add_argument("--gtdb-bac-meta", type=Path, default=None)
    ap.add_argument("--gtdb-arc-meta", type=Path, default=None)
    ap.add_argument("--gtdb-genomes-tar", type=Path, default=None)
    ap.add_argument("--kmer-len", type=int, default=35)
    ap.add_argument("--minimizer-len", type=int, default=31)
    ap.add_argument("--minimizer-spaces", type=int, default=7)
    ap.add_argument("--repo-root", type=Path, default=Path.cwd())
    args = ap.parse_args()

    hash_k2d = args.db_dir / "hash.k2d"
    if not hash_k2d.exists():
        sys.exit(f"Error: {hash_k2d} does not exist; refusing to write manifest for an absent index")

    accessions = []
    for raw in args.accessions_file.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        accessions.append(line.split()[0])

    manifest = {
        "name": args.name,
        "label": args.label,
        "built_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "user": os.environ.get("USER", "unknown"),
        "git_commit": git_commit(args.repo_root),
        "gtdb_release": args.gtdb_release or None,
        "accessions_file": file_record(args.accessions_file, "accessions"),
        "accessions": accessions,
        "supplemental_taxonomy": file_record(args.supplemental_taxonomy, "supplemental_taxonomy"),
        "supplemental_urls": file_record(args.supplemental_urls, "supplemental_urls"),
        "merged_taxonomy": file_record(args.merged_taxonomy, "merged_taxonomy"),
        "gtdb_sources": {
            "bac120_taxonomy": file_record(args.gtdb_bac_tax, "bac120_taxonomy"),
            "ar53_taxonomy": file_record(args.gtdb_arc_tax, "ar53_taxonomy"),
            "bac120_metadata": file_record(args.gtdb_bac_meta, "bac120_metadata"),
            "ar53_metadata": file_record(args.gtdb_arc_meta, "ar53_metadata"),
            "genomes_reps_tar": file_record(args.gtdb_genomes_tar, "genomes_reps_tar"),
        },
        "kraken2_index": {
            "hash_k2d": file_record(hash_k2d, "hash_k2d"),
            "opts_k2d": file_record(args.db_dir / "opts.k2d", "opts_k2d"),
            "taxo_k2d": file_record(args.db_dir / "taxo.k2d", "taxo_k2d"),
        },
        "tools": {
            "kraken2": get_version([shutil.which("kraken2") or "kraken2", "--version"]),
            "kraken2-build": get_version([shutil.which("kraken2-build") or "kraken2-build", "--version"]),
            "ncbi-datasets-cli": get_version([shutil.which("datasets") or "datasets", "--version"]),
            "taxonkit": get_version([shutil.which("taxonkit") or "taxonkit", "version"]),
        },
        "build_params": {
            "kmer_len": args.kmer_len,
            "minimizer_len": args.minimizer_len,
            "minimizer_spaces": args.minimizer_spaces,
        },
    }

    out = args.db_dir / f"MANIFEST.{args.name}.json"
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
