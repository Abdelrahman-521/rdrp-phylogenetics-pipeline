// Maximum-likelihood trees and the side-by-side comparison.

process IQTREE {
    tag "${tool}"
    label 'tree'
    publishDir "${params.outdir}/04_trees", mode: 'copy', pattern: "${tool}.*"

    input:
    tuple val(tool), path(alignment)

    output:
    tuple val(tool), path(alignment), path("${tool}.treefile"), path("${tool}.iqtree"), emit: result
    path "${tool}.log", emit: log

    script:
    """
    iqtree2 -s ${alignment} --prefix ${tool} \\
        -m ${params.iqtree_model} ${params.iqtree_extra} \\
        -bb ${params.bootstraps} -nt ${task.cpus} -seed ${params.seed} -redo -quiet
    """
}

process COMPARE_TREES {
    tag "ClustalW vs T-Coffee"
    publishDir "${params.outdir}/comparison", mode: 'copy'

    input:
    path clustal_aln
    path tcoffee_aln
    path clustal_tree
    path tcoffee_tree
    path clustal_report
    path tcoffee_report

    output:
    path 'comparison.md'
    path 'comparison.tsv'
    path 'comparison_tanglegram.png'

    script:
    """
    compare_trees.py \\
        --clustal-aln ${clustal_aln} --tcoffee-aln ${tcoffee_aln} \\
        --clustal-tree ${clustal_tree} --tcoffee-tree ${tcoffee_tree} \\
        --clustal-report ${clustal_report} --tcoffee-report ${tcoffee_report} \\
        --prefix comparison
    """
}
