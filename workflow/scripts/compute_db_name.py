#!/usr/bin/env python3
"""Compute deterministic kraken2 database directory name from a release + label + accessions file.

Usage:
    compute_db_name.py --accessions <path> --label <slug> [--gtdb-release <ver>]
    compute_db_name.py --accessions <path> --label <slug> [--gtdb-release <ver>] --print-accessions
"""
from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path


_ACC_RE = re.compile(r"^(GC[AF])_\d{9}\.\d+$")


def read_accessions(path: Path) -> list[str]:
    accs: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        token = line.split()[0]
        if not _ACC_RE.match(token):
            sys.exit(
                f"Error: {path}: '{token}' is not a valid versioned NCBI accession "
                "(expected GCF_NNNNNNNNN.V or GCA_NNNNNNNNN.V)"
            )
        accs.append(token)
    if not accs:
        sys.exit(f"Error: {path}: no accessions found")
    return accs


def canonical(accessions: list[str]) -> list[str]:
    return sorted(set(accessions))


def short_sha(accessions: list[str]) -> str:
    blob = ("\n".join(canonical(accessions)) + "\n").encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:8]


def slugify(label: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "-", label).strip("-")
    if not s:
        sys.exit(f"Error: label '{label}' produced an empty slug")
    return s


def compute_name(label: str, accessions: list[str], gtdb_release: str | None) -> str:
    sha = short_sha(accessions)
    slug = slugify(label)
    if gtdb_release:
        return f"kraken2-gtdb-r{gtdb_release}-{slug}-{sha}"
    return f"kraken2-{slug}-{sha}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", type=Path, required=True)
    ap.add_argument("--label", type=str, required=True)
    ap.add_argument("--gtdb-release", type=str, default=None)
    ap.add_argument("--print-accessions", action="store_true",
                    help="Print canonical accession list (one per line) to stdout instead of the name")
    args = ap.parse_args()

    accessions = read_accessions(args.accessions)
    if args.print_accessions:
        print("\n".join(canonical(accessions)))
    else:
        print(compute_name(args.label, accessions, args.gtdb_release))
    return 0


if __name__ == "__main__":
    sys.exit(main())
