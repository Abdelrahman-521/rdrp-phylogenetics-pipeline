#!/bin/bash

# This is a script used to download the accession files for use in the subsequent pipeline.
# The accessions are saved to a file called sequences.fasta (in fasta format)

# Define the list of accession numbers
accessions=(
    "MG833031.1" "CAA25029.1" "BAA03151.1" "AAP97137.1" "NP_041277.1" "CAA25416.1" "NP_056777.1"
    "BAA00168.1" "ABS84820.1" "YP_001285409.1" "NP_853560.2" "YP_145791.1" "AAQ64627.1"
    "BAD06930.1" "NP_733845.1" "NP_277061.1" "NP_620564.1" "AAF80998.1" "AAN63804.2"
    "AAF00472.1" "BAA32553.1" "NP_049374.1" "BAA21898" "AAC58807.1" "NP_046155.1"
    "ADP24157.1" "YP_009337666.1" "YP_009337152.1" "YP_009336629.1" "ASK12212.1"
    "APG77945.1" "ASK12194.1" "ASK12217.1" "APG77984.1" "LC310707"
    "KU645789" "ABB17263.2" "YP_004063985.1"
)

# Output file for sequences in FASTA format
output_file="sequences.fasta"

# Loop through each accession number and download the sequence
for acc in "${accessions[@]}"; do
    echo "Fetching $acc..."
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$acc&rettype=fasta&retmode=text" >> "$output_file"
    sleep 1  # avoid rate-limiting
done

echo "All sequences downloaded and saved to $output_file."