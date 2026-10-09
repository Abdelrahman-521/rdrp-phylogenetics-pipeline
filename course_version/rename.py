# Define a dictionary mapping full headers to abbreviations
header_to_abbreviation = {
    "CAA25029.1_1178_1670 SSF56672": "CPMV",
    "BAA03151.1_2260_2777 SSF56672": "PYFV",
    "AAP97137.1_1221_1664 SSF56672": "HaRNAV",
    "NP_041277.1_1750_2208 SSF56672": "PV",
    "CAA25416.1_1869_2331 SSF56672": "FMDV",
    "NP_056777.1_1838_2291 SSF56672": "EMCV",
    "BAA00168.1_1701_2156 SSF56672": "HRV",
    "ABS84820.1_2447_2956 SSF56672": "SBPV",
    "YP_001285409.1_2470_2969 SSF56672": "BrBV",
    "NP_853560.2_2380_2889 SSF56672": "DWV",
    "YP_145791.1_2380_2889 SSF56672": "VDV-1",
    "AAQ64627.1_2452_2931 SSF56672": "EoV",
    "BAD06930.1_2380_2889 SSF56672": "KV",
    "NP_733845.1_1491_1997 SSF56672": "ALPV",
    "NP_277061.1_2451_2929 SSF56672": "PnV",
    "NP_620564.1_1139_1626 SSF56672": "BQCV",
    "AAF80998.1_1241_1765 SSF56672": "CrPV",
    "AAN63804.2_1391_1899 SSF56672": "ABPV",
    "AAF00472.1_1267_1778 SSF56672": "TrV",
    "BAA32553.1_1271_1770 SSF56672": "HiPV",
    "NP_049374.1_2339_2847 SSF56672": "SBV",
    "BAA21898.1_1297_1795 SSF56672": "PSIV",
    "AAC58807.1_1230_1748 SSF56672": "DCV",
    "NP_046155.1_1449_1965 SSF56672": "RhPV",
    "ADP24157.1_2455_2949 SSF56672": "IFV",
    "YP_009337666.1_2332_2896 SSF56672": "HplV-35",
    "YP_009337152.1_2311_2848 SSF56672": "HplV-34",
    "YP_009336629.1_2300_2816 SSF56672": "HuAV 1",
    "ASK12212.1_1647_2131 SSF56672": "LniV-1",
    "APG77945.1_1800_2297 SSF56672": "ShiV-8",
    "ASK12194.1_1780_2279 SSF56672": "SINV-4",
    "ASK12217.1_1635_2119 SSF56672": "SINV-2",
    "APG77984.1_1681_2186 SSF56672": "HplV-81",
    "ABB17263.2_1599_2081 SSF56672": "TSV",
    "YP_004063985.1_1669_2145 SSF56672": "MCDV",
    "lcl|MG833031.1_prot_AWC26954.1_1_2530_3088 SSF56672": "CuPV-1",
    "lcl|LC310707.1_prot_BBC14926.1_1_2237_2725 SSF56672": "AaIFV",
    "lcl|KU645789.1_prot_AOT85373.1_1_2535_3035 SSF56672": "MV",
}

# Input and output file paths
input_fasta = "SSF56672_sequences.fasta"
output_fasta = "updated_headers_full_match.fasta"

# Process the FASTA file
with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            # Extract the full header and map it to the abbreviation
            full_header = line[1:].strip()  # Remove ">" and newline
            abbreviation = header_to_abbreviation.get(full_header, full_header)
            new_header = f">{abbreviation}\n"  # Write the updated header
            outfile.write(new_header)
        else:
            outfile.write(line)  # Write sequence lines unchanged

print(f"Updated FASTA file written to: {output_fasta}")
