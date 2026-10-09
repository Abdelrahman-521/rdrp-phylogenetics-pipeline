// Getting from a list of accessions to clean, named RdRp domain sequences.

process FETCH_PROTEINS {
    tag "NCBI protein"
    publishDir "${params.outdir}/01_download", mode: 'copy'
    // NCBI drops connections now and then: wait and try again.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path accessions

    output:
    path 'sequences.fasta', emit: fasta
    path 'fetch_log.tsv',   emit: log

    script:
    """
    fetch_proteins.py ${accessions} -o sequences.fasta --log fetch_log.tsv
    """
}

process REPAIR_FASTA {
    tag "clean + recover"
    publishDir "${params.outdir}/01_download", mode: 'copy'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path raw
    path accessions

    output:
    path 'fixed_sequences.fasta'

    script:
    """
    repair_fasta.py ${raw} ${accessions} -o fixed_sequences.fasta
    """
}

process EXTRACT_DOMAIN {
    tag "SSF56672"
    publishDir "${params.outdir}/02_domains", mode: 'copy'

    input:
    path fasta
    path regions

    output:
    path 'rdrp_domains.fasta'

    script:
    """
    extract_domain.py ${fasta} ${regions} -o rdrp_domains.fasta --format ${params.regions_format}
    """
}

process RENAME_HEADERS {
    tag "short names"
    publishDir "${params.outdir}/02_domains", mode: 'copy'

    input:
    path domains
    path names

    output:
    path 'named_domains.fasta'

    script:
    """
    rename_headers.py ${domains} ${names} -o named_domains.fasta
    """
}
