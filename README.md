# Viral RdRp Phylogenetics Pipeline

A reproducible pipeline that builds evolutionary trees for **RNA-dependent RNA polymerase (RdRp)**, the enzyme that RNA viruses use to copy their genomes. It downloads 38 viral protein sequences from NCBI, repairs the download, aligns the sequences two different ways, and compares the maximum-likelihood trees that come out of each alignment.

Group project for COMP 3550 / BIOL 3951 (Bioinformatics) at Memorial University, Fall 2024.

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1. Download | `scripts/getseqs.sh` | Fetches 38 accessions from NCBI E-utilities (`efetch`), with a pause between calls to respect rate limits |
| 2. Repair | `scripts/checkerrors.py` | Strips error output from the FASTA file and re-downloads any accession that failed. Some IDs live in the nucleotide database rather than the protein database, so it retries against `nuccore` |
| 3. Clean | `scripts/cleanrdrp.py`, `scripts/trimseqs.py` | Normalises headers and trims sequences for alignment |
| 4. Align | `scripts/clustal_rdrp.sh`, `scripts/tcoffee_rdrp.sh` | Aligns with **ClustalW** and **T-Coffee** in parallel, so each alignment can cross-check the other |
| 5. Domain filter | `scripts/get_superfamily.py`, `scripts/rename.py` | Keeps the sequences that InterPro assigns to the RdRp superfamily (SSF56672) and gives them readable names |
| 6. Trees | IQ-TREE 2 | Picks the best substitution model automatically (`-m MFP`) and runs 1,000 ultrafast bootstraps (`-bb 1000`) |
| 7. Orchestration | `scripts/main.nf` | Nextflow entry point that chains the steps |

Exact commands and tool versions are in [`docs/pipeline_order_and_versions.txt`](docs/pipeline_order_and_versions.txt).

## Results

| ClustalW tree | T-Coffee tree |
|---|---|
| ![ClustalW tree](data/clustaltree.png) | ![T-Coffee tree](data/tcoffeetree.png) |

- Alignments: `results/rdrp_clustal_aligned.fasta`, `results/rdrp_tcoffee_aligned.fasta`
- IQ-TREE reports and Newick trees: `results/iqtree/`
- Slides: [final presentation](docs/project_presentation.pdf) and [proposal](docs/proposal_presentation.pdf)

## Run it

```bash
conda install -c bioconda clustalw t-coffee iqtree biopython nextflow
bash scripts/getseqs.sh            # writes sequences.fasta
python scripts/checkerrors.py      # repairs failed downloads
bash scripts/clustal_rdrp.sh && bash scripts/tcoffee_rdrp.sh
iqtree2 -s clustalALN.fasta -m MFP -bb 1000 -nt AUTO
```

## Tech
Bash · Python · Biopython · NCBI E-utilities · ClustalW · T-Coffee · IQ-TREE 2 · InterPro · Nextflow · FigTree

## Credits
Built as a team project. My steps are logged in [`docs/pipeline_order_and_versions.txt`](docs/pipeline_order_and_versions.txt): the InterPro superfamily filter, renaming, the ClustalW and T-Coffee re-alignments, and the IQ-TREE runs. Thanks to my teammates and the course instructors.
