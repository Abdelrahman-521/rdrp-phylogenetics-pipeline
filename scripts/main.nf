nextflow.enable.dsl=2

process hello {
    input:
    val name

    output:
    path 'greeting.txt'

    script:
    """
    echo "Hello, $name!" > greeting.txt
    """
}

workflow {
    hello("World")
}
