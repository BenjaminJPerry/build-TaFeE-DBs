#!/usr/bin/env python3
"""Resolve a GTDB release identifier (e.g. '220.0') into the 5 canonical download URLs."""
from __future__ import annotations


GTDB_BASE = "https://data.ace.uq.edu.au/public/gtdb/data/releases"


def gtdb_urls(release: str) -> dict[str, str]:
    """Return the 5 canonical GTDB URLs for the given release (e.g. '220.0')."""
    major = release.split(".", 1)[0]
    base = f"{GTDB_BASE}/release{major}/{release}"
    return {
        "bac_tax_url":  f"{base}/bac120_taxonomy_r{major}.tsv.gz",
        "arc_tax_url":  f"{base}/ar53_taxonomy_r{major}.tsv.gz",
        "bac_meta_url": f"{base}/bac120_metadata_r{major}.tsv.gz",
        "arc_meta_url": f"{base}/ar53_metadata_r{major}.tsv.gz",
        "genomes_url":  f"{base}/genomic_files_reps/gtdb_genomes_reps_r{major}.tar.gz",
    }


if __name__ == "__main__":
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("release")
    args = ap.parse_args()
    print(json.dumps(gtdb_urls(args.release), indent=2))
