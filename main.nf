#!/usr/bin/env nextflow
/*
 * PhyloBot: does the alignment tool change the evolutionary tree?
 *
 * Downloads viral RdRp proteins from NCBI, repairs the download, cuts out the
 * RdRp domain, aligns it with ClustalW and T-Coffee side by side, builds a
 * maximum-likelihood tree from each alignment with IQ-TREE, then measures how
 * different the two trees are.
 *
 *   nextflow run . -profile local          # full run from NCBI
 *   nextflow run . -profile local,test     # offline run from bundled domain sequences
 */

nextflow.enable.dsl = 2

include { FETCH_PROTEINS; REPAIR_FASTA; EXTRACT_DOMAIN; RENAME_HEADERS } from './modules/prepare'
include { ALIGN_CLUSTALW; ALIGN_TCOFFEE }                                from './modules/align'
include { IQTREE; COMPARE_TREES }                                        from './modules/trees'

def summary() {
    log.info """
    P H Y L O B O T
    ===============
    accessions : ${params.domains ? '(skipped: starting from --domains)' : params.accessions}
    domains    : ${params.domains ?: 'cut from NCBI downloads using ' + params.regions}
    names      : ${params.names}
    bootstraps : ${params.bootstraps}   seed: ${params.seed}   model: ${params.iqtree_model} ${params.iqtree_extra}
    outdir     : ${params.outdir}
    """.stripIndent()
}

workflow {
    summary()

    names = file(params.names, checkIfExists: true)

    if (params.domains) {
        domains = Channel.fromPath(params.domains, checkIfExists: true)
    } else {
        accessions = file(params.accessions, checkIfExists: true)
        raw        = FETCH_PROTEINS(accessions)
        fixed      = REPAIR_FASTA(raw.fasta, accessions)
        domains    = EXTRACT_DOMAIN(fixed, file(params.regions, checkIfExists: true))
    }

    named = RENAME_HEADERS(domains, names)

    // The two aligners run in parallel on the same input.
    alignments = ALIGN_CLUSTALW(named).mix(ALIGN_TCOFFEE(named))

    trees = IQTREE(alignments).result

    // Collect both trees, then compare them in one step.
    by_tool = trees.map { tool, aln, tree, report -> [tool, [aln, tree, report]] }
                   .collect(flat: false)
                   .map { pairs -> pairs.collectEntries() }

    COMPARE_TREES(
        by_tool.map { it['clustalw'][0] }, by_tool.map { it['tcoffee'][0] },
        by_tool.map { it['clustalw'][1] }, by_tool.map { it['tcoffee'][1] },
        by_tool.map { it['clustalw'][2] }, by_tool.map { it['tcoffee'][2] }
    )
}

workflow.onComplete {
    log.info(workflow.success
        ? "\nDone in ${workflow.duration}. Report: ${params.outdir}/comparison/comparison.md\n"
        : "\nFailed: ${workflow.errorMessage}\n")
}
