# Course version (Fall 2024)

The scripts exactly as they were at the end of COMP 3550 / BIOL 3951. The team ran them by hand, one after another:

| Order | File | What it did |
|---|---|---|
| 1 | `getseqs.sh` | Downloaded the 38 proteins from NCBI, one `curl` call each |
| 2 | `checkerrors.py` | Removed NCBI error text from the download and re-fetched the three accessions that are nucleotide records, as translated CDS from `nuccore` |
| 3 | InterProScan (web) | Annotated the proteins and produced `interpro_results.tsv` (not kept) |
| 4 | `get_superfamily.py` | Cut out the SSF56672 (RdRp) domain from each protein |
| 5 | `rename.py` | Replaced accessions with short virus names (`named_seq.fasta`) |
| 6 | ClustalW and T-Coffee | Aligned the domains (commands in `pipeline_order_and_versions.txt`) |
| 7 | IQ-TREE 2.3.6 | One maximum-likelihood tree per alignment (`results/iqtree/`) |
| 8 | FigTree | Compared the two trees by eye |

`trimseqs.py`, `cleanrdrp.py`, `clustal.sh` and the two `results/rdrp_*_aligned.fasta` files are from the first plan: re-doing the CuPV-1 paper's MEGA analysis on a consensus-trimmed region. Our ClustalW tree didn't match the paper's, so the team moved to the InterPro domain approach above.

`main.nf` and `nextflow.config` were a Nextflow "Hello World" test. Wiring the steps into Nextflow was listed as future work on the final slide; it's done now in the top-level `main.nf`.

I (Abdelrahman Conber) wrote the Python scripts here. Team: Abdelrahman Conber, Claire Gallant and Haley Leonard.
