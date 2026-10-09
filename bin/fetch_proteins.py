#!/usr/bin/env python3
"""Step 1: download protein sequences from NCBI.

Pipeline version of the course's getseqs.sh. Asks NCBI's protein database for
every accession, one request at a time with a pause between calls. Accessions
the protein database doesn't have are logged, not fatal: the repair step
(repair_fasta.py) fetches them from the nucleotide database instead.
"""

import argparse
import os
import sys

from phylobot_utils import efetch, parse_fasta_text, read_list, write_fasta


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accessions", help="text file, one accession per line")
    parser.add_argument("-o", "--output", default="sequences.fasta")
    parser.add_argument("--log", default="fetch_log.tsv")
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    args = parser.parse_args()

    api_key = os.environ.get("NCBI_API_KEY")
    accessions = read_list(args.accessions)
    records = []
    with open(args.log, "w") as log:
        log.write("accession\tdatabase\tstatus\trecords\n")
        for acc in accessions:
            text = efetch("protein", acc, "fasta", email=args.email, api_key=api_key)
            got, junk = parse_fasta_text(text)
            status = "ok" if got else "not_in_protein_db"
            print(f"{acc}: {status}" + (f" ({len(junk)} junk lines dropped)" if junk else ""), file=sys.stderr)
            log.write(f"{acc}\tprotein\t{status}\t{len(got)}\n")
            records.extend(got)

    write_fasta(records, args.output)
    print(f"Saved {len(records)} of {len(accessions)} sequences to {args.output}", file=sys.stderr)
    return 0 if records else 1


if __name__ == "__main__":
    sys.exit(main())
