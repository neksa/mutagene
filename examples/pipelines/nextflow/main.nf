nextflow.enable.dsl=2
params.outdir = 'results'

process DEMO {
    output:
    path 'demo/matrix.tsv', emit: matrix
    script:
    '''
    mutagene demo --output-dir demo
    '''
}

process SIGNATURE {
    publishDir params.outdir, mode: 'copy'
    input:
    path matrix
    output:
    path 'signatures.tsv'
    path 'diagnostics.json'
    script:
    """
    mutagene signature -f MATRIX -i '${matrix}' -s COSMICv2 \\
        -o signatures.tsv --diagnostics-out diagnostics.json
    """
}

workflow {
    DEMO()
    SIGNATURE(DEMO.out.matrix)
}
