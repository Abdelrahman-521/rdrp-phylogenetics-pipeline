from Bio import AlignIO

# input / output file
alignment_file = "aligned_sequences.fasta" 
output_file = "trimmed_rdrp.fasta"

# Define the consensus-based start and end RdRp position (1-based indexing)
consensus_start = 3296
consensus_end = 3648

# Read the alignment
alignment = AlignIO.read(alignment_file, "fasta")

# Open output file for writing
with open(output_file, "w") as outfile:
    for record in alignment:
        # Extract the region corresponding to the consensus positions
        trimmed_seq = record.seq[consensus_start - 1:consensus_end]
        # Write the trimmed sequence to the output file
        outfile.write(f">{record.id}\n{trimmed_seq}\n")

