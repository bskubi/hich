include {source; emptyOnLastStep} from './extraops.nf'

process ChromsizesProc {
    publishDir params.general.publish.chromsizes ? params.general.publish.chromsizes : "results",
               saveAs: {params.general.publish.chromsizes ? it : null},
               mode: params.general.publish.mode

    conda "bioconda::ucsc-fasize"
    container 'quay.io/biocontainers/ucsc-fasize:332--0'
    label 'smallResource'
    memory 8.GB

    input:
    tuple path(reference), val(assembly), val(chromsizes)

    output:
    tuple val(assembly), path(chromsizes)

    shell:
    "faSize -detailed -tab ${reference} > ${chromsizes}"

    stub:
    "touch ${chromsizes}"
}

workflow Chromsizes {
    take:
    samples

    main:

    source(ChromsizesProc,
           samples,
           "chromsizes",
           ["genomeReference", "assembly", "chromsizes"],
           ["assembly", "chromsizes"],
           {"${it.assembly}.sizes"},
           "assembly",
            {true}) | set{samples}
    
    samples = emptyOnLastStep("Chromsizes", samples)

    emit:
    samples
}