// The comparison at the heart of the project: same sequences, two aligners.

process ALIGN_CLUSTALW {
    tag "ClustalW"
    label 'align'
    publishDir "${params.outdir}/03_alignments", mode: 'copy'

    input:
    path fasta

    output:
    tuple val('clustalw'), path('clustalw_aligned.fasta')

    script:
    """
    clustalw -INFILE=${fasta} -OUTFILE=clustalw_aligned.fasta -OUTPUT=FASTA -OUTORDER=INPUT > clustalw.log
    """
}

process ALIGN_TCOFFEE {
    tag "T-Coffee"
    label 'align'
    publishDir "${params.outdir}/03_alignments", mode: 'copy'

    input:
    path fasta

    output:
    tuple val('tcoffee'), path('tcoffee_aligned.fasta')

    script:
    """
    # T-Coffee writes caches and temp files under \$HOME; keep them in the task folder.
    export HOME_4_TCOFFEE=\$PWD/.t_coffee TMP_4_TCOFFEE=\$PWD/.t_coffee/tmp CACHE_4_TCOFFEE=\$PWD/.t_coffee/cache
    mkdir -p \$TMP_4_TCOFFEE \$CACHE_4_TCOFFEE
    t_coffee -infile ${fasta} -output fasta_aln -outfile tcoffee_aligned.fasta \\
        -n_core ${task.cpus} -quiet > tcoffee.log 2>&1
    """
}
