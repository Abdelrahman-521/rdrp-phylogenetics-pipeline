#!/usr/bin/env python3
"""Step 4: give every domain a short virus name so the trees are readable.

Pipeline version of the course's rename.py. The mapping lives in
data/virus_names.tsv instead of inside the code. Spaces become underscores
because tree tools cut names at the first space (the course run turned
'HuAV 1' into 'HuAV').
"""

import argparse
import re
import sys

from phylobot_utils import parse_fasta_text, read_tsv, record_id, write_fasta


def short_name_for(rec_id, names):
    seq_id = re.sub(r"_\d+_\d+$", "", rec_id)
    return names.get(seq_id)


def rename(records, names):
    out, unnamed, used = [], [], set()
    for header, seq in records:
        rid = record_id(header)
        name = short_name_for(rid, names)
        if name is None:
            unnamed.append(rid)
            name = rid
        name = re.sub(r"\s+", "_", name.strip())
        if name in used:
            raise ValueError(f"two sequences would both be called {name!r}")
        used.add(name)
        out.append((name, seq))
    return out, unnamed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta")
    parser.add_argument("names", help="TSV with seq_id and short_name columns")
    parser.add_argument("-o", "--output", default="named_domains.fasta")
    args = parser.parse_args()

    names = {row["seq_id"]: row["short_name"] for row in read_tsv(args.names)}
    with open(args.fasta) as handle:
        records, _ = parse_fasta_text(handle.read())
    renamed, unnamed = rename(records, names)
    write_fasta(renamed, args.output)
    for rid in unnamed:
        print(f"WARNING no short name for {rid}; kept the accession", file=sys.stderr)
    print(f"Renamed {len(renamed)} sequences to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
