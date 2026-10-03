from Bio import SeqIO

# input / output
input_file = "trimmed_rdrp.fasta"
output_file = "rdrp_regions_cleaned.fasta"

# open output file for writing, write trimmed file to new file replacing gaps (-) with ""
with open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        # Remove gaps from the sequence using string replacement
        clean_sequence = str(record.seq).replace("-", "")
        # Write the cleaned sequence to the new file
        output_handle.write(f">{record.id}\n{clean_sequence}\n")

