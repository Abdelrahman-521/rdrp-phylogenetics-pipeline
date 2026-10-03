# This is the error checking script.
# It cleans up the sequence file, by removing download error output and redownloads missed sequences (nt sequences)
# This accounts for some of the original accessions being on the nuccore database instead of the protein database.

#In order to download missed sequences using bash's curl (within python) we need to import package os.
import os

# Paths to the original FASTA file and the new final output file
original_fasta_file = "sequences.fasta"
final_output_file = "fixed_sequences.fasta"

# Full list of accession numbers
accessions = [
    "MG833031.1", "CAA25029.1", "BAA03151.1", "AAP97137.1", "NP_041277.1", "CAA25416.1", "NP_056777.1",
    "BAA00168.1", "ABS84820.1", "YP_001285409.1", "NP_853560.2", "YP_145791.1", "AAQ64627.1",
    "BAD06930.1", "NP_733845.1", "NP_277061.1", "NP_620564.1", "AAF80998.1", "AAN63804.2",
    "AAF00472.1", "BAA32553.1", "NP_049374.1", "BAA21898", "AAC58807.1", "NP_046155.1",
    "ADP24157.1", "YP_009337666.1", "YP_009337152.1", "YP_009336629.1", "ASK12212.1",
    "APG77945.1", "ASK12194.1", "ASK12217.1", "APG77984.1", "LC310707",
    "KU645789", "ABB17263.2", "YP_004063985.1"
]

# Function 1: cleans original file of error outputs, which start with "+"
def clean_fasta_file(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    cleaned_lines = [line for line in lines if not line.startswith("+")]
    with open(file_path, "w") as f:
        f.writelines(cleaned_lines)

# Function 2: checks if all accessions from the variable "accessions" are in the file. The missing ones are returned by the function.
def check_accessions(file_path, accessions):
    with open(file_path, "r") as f:
        content = f.read()
    missing_accessions = [acc for acc in accessions if f">{acc}" not in content]
    return missing_accessions

# Function 3: downloads missing accessions to end of file using the nuccore database.
# We are assuming that if getseqs.sh missed a download, it is due to the accession not being in the protein database and is instead in the nuccore database.
def download_missing_sequences(missing_accessions, final_file):
    for accession in missing_accessions:
        print(f"Downloading protein translation for CDS of {accession}...")
        response = os.system(
            f'curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta_cds_aa&retmode=text" >> "{final_file}"'
        )
        if response == 0:
            print(f"Completed download for {accession}.")
        else:
            print(f"Error downloading {accession}")

# Main function
def main():
    # Function 1: Clean
    clean_fasta_file(original_fasta_file)

    # Function 2: Make a list of any missing accessions
    missing_accessions = check_accessions(original_fasta_file, accessions)

    # Create a cleaned "final" file (not added missing accessions yet)
    with open(final_output_file, "w") as final:
        with open(original_fasta_file, "r") as original:
            final.writelines(original.readlines())

    # Function 3: Download the missing accessions from the nuccore database to append to the final file
    if missing_accessions:
        download_missing_sequences(missing_accessions, final_output_file)
        print(f"Missing sequences added to '{final_output_file}'.")

    print(f"Error-checking completed. Final sequences saved to '{final_output_file}'.")

# Run the main function if the script is called. 
if __name__ == "__main__":
    main()
