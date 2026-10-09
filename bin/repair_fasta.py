#!/usr/bin/env python3
"""Step 2: clean the download and recover anything that's missing.

Pipeline version of the course's checkerrors.py, which Abdel wrote:
  1. drop anything in the file that isn't a real FASTA record (NCBI error text),
  2. find which accessions are still missing,
  3. fetch those from the nucleotide database (nuccore) as translated CDS,
     because some of the paper's accessions are nucleotide records,
  4. fail loudly if anything is still missing, so the workflow can retry
     instead of quietly building a tree from fewer sequences.
"""

import argparse
import os
import sys

from phylobot_utils import (
    efetch,
    missing_accessions,
    parse_fasta_text,
    read_list,
    write_fasta,
)


def repair(text, accessions, fetch):
    """Return (records, report). `fetch(accession)` returns raw FASTA text."""
    records, junk = parse_fasta_text(text)
    report = {"junk_lines": len(junk), "recovered": [], "still_missing": []}
    for acc in missing_accessions(records, accessions):
        got, _ = parse_fasta_text(fetch(acc))
        if got:
            records.extend(got)
            report["recovered"].append(acc)
    report["still_missing"] = missing_accessions(records, accessions)
    return records, report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="raw download from fetch_proteins.py")
    parser.add_argument("accessions", help="text file, one accession per line")
    parser.add_argument("-o", "--output", default="fixed_sequences.fasta")
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    args = parser.parse_args()

    api_key = os.environ.get("NCBI_API_KEY")

    def fetch_cds(acc):
        print(f"{acc}: fetching translated CDS from nuccore", file=sys.stderr)
        return efetch("nuccore", acc, "fasta_cds_aa", email=args.email, api_key=api_key)

    with open(args.fasta) as handle:
        records, report = repair(handle.read(), read_list(args.accessions), fetch_cds)

    write_fasta(records, args.output)
    print(
        f"Dropped {report['junk_lines']} junk lines; recovered {len(report['recovered'])} "
        f"accessions from nuccore ({', '.join(report['recovered']) or 'none'}); "
        f"{len(records)} records saved to {args.output}",
        file=sys.stderr,
    )
    if report["still_missing"]:
        print("Still missing: " + ", ".join(report["still_missing"]), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
