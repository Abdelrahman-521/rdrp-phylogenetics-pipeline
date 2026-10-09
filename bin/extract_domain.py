#!/usr/bin/env python3
"""Step 3: cut out the RdRp domain from each full-length protein.

Pipeline version of the course's get_superfamily.py. Takes domain coordinates
either from this repo's simple regions table or straight from an InterProScan
TSV, keeps the SUPERFAMILY SSF56672 hits (the RNA-dependent RNA polymerase
fold) and writes one domain sequence per protein.
"""

import argparse
import sys
from collections import defaultdict

from phylobot_utils import parse_fasta_text, read_tsv, record_id, write_fasta


def load_regions(path, fmt, signature):
    regions = defaultdict(list)
    if fmt == "interproscan":
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 7 and fields[4] == signature:
                    regions[fields[0]].append((int(fields[6]), int(fields[7])))
    else:
        for row in read_tsv(path):
            if row["signature"] == signature:
                regions[row["seq_id"]].append((int(row["start"]), int(row["end"])))
    return regions


def extract(records, regions, signature):
    out, problems = [], []
    seen = set()
    for header, seq in records:
        rid = record_id(header)
        for start, end in regions.get(rid, []):
            if end > len(seq):
                problems.append(f"{rid}: region {start}-{end} is past the end ({len(seq)} aa)")
                continue
            out.append((f"{rid}_{start}_{end} {signature}", seq[start - 1 : end]))
            seen.add(rid)
    for rid in regions:
        if rid not in seen and not any(rid in p for p in problems):
            problems.append(f"{rid}: listed in the regions file but not found in the FASTA")
    return out, problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta")
    parser.add_argument("regions")
    parser.add_argument("-o", "--output", default="rdrp_domains.fasta")
    parser.add_argument("--format", choices=["regions", "interproscan"], default="regions")
    parser.add_argument("--signature", default="SSF56672")
    args = parser.parse_args()

    with open(args.fasta) as handle:
        records, _ = parse_fasta_text(handle.read())
    domains, problems = extract(records, load_regions(args.regions, args.format, args.signature), args.signature)
    write_fasta(domains, args.output)
    for problem in problems:
        print("WARNING " + problem, file=sys.stderr)
    print(f"Extracted {len(domains)} {args.signature} domains to {args.output}", file=sys.stderr)
    return 0 if domains and not problems else 1


if __name__ == "__main__":
    sys.exit(main())
