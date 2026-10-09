from Bio import SeqIO

# File paths
fasta_file = "fixed_sequences.fasta"  
interpro_file = "interpro_results.tsv"  
output_file = "SSF56672_sequences.fasta"

# Define the superfamily to extract
target_superfamily = "SSF56672"

# Parse InterPro results to find relevant regions
regions_to_extract = {}
with open(interpro_file, "r") as interpro:
    for line in interpro:
        fields = line.strip().split("\t")
        seq_id, superfamily, start, end = fields[0], fields[4], int(fields[6]), int(fields[7])
        if superfamily == target_superfamily:
            if seq_id not in regions_to_extract:
                regions_to_extract[seq_id] = []
            regions_to_extract[seq_id].append((start, end))

# Extract the sequences from the FASTA file
with open(output_file, "w") as output:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in regions_to_extract:
            for start, end in regions_to_extract[record.id]:
                # Extract the rdrp region
                extracted_seq = record.seq[start-1:end]
                output.write(f">{record.id}_{start}_{end} SSF56672\n{extracted_seq}\n")

print(f"Extracted SSF56672 sequences saved to {output_file}.")
