include {source; emptyOnLastStep; pack2; isExistingFile} from './extraops.nf'

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
    
    samples
        | filter {!isExistingFile(it.chromsizes)}
        | map{tuple(it.genomeReference, it.assembly, "${it.assembly}.sizes")}
        | unique
        | ChromsizesProc
        | map{assembly, chromsizes -> [assembly: assembly, chromsizes: chromsizes]}
        | set{result}
    pack2(samples, result, "assembly") | set{samples}

    samples = emptyOnLastStep("Chromsizes", samples)

    emit:
    samples
}