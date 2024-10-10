include {emptyOnLastStep; pack; skip; isExistingFile} from './extraops.nf'

process ChromsizesProc {
    publishDir params.general.publish.chromsizes ? params.general.publish.chromsizes : "results",
               saveAs: {params.general.publish.chromsizes ? it : null},
               mode: params.general.publish.mode

    conda "bioconda::ucsc-fasize"
    container 'quay.io/biocontainers/ucsc-fasize:332--0'
    label 'smallResource'
    memory 8.GB

    input:
    tuple val(genomeReferenceString), path(genomeReference), val(assembly), val(chromsizes)

    output:
    tuple val(genomeReferenceString), val(assembly), path(chromsizes)

    shell:
    "faSize -detailed -tab ${genomeReference} > ${chromsizes}"

    stub:
    "touch ${chromsizes}"
}

workflow Chromsizes {
    take:
    samples

    main:
    
    samples
        | filter {!skip("chromsizes") && !isExistingFile(it.chromsizes)}
        | map{tuple(it.genomeReference, it.genomeReference, it.assembly, "${it.assembly}.sizes")}
        | unique
        | ChromsizesProc
        | map{genomeReference, assembly, chromsizes -> [genomeReference: file(genomeReference), assembly: assembly, chromsizes: chromsizes]}
        | set{result}
    pack(samples, result, ["genomeReference", "assembly"]) | set{samples}

    samples = emptyOnLastStep("chromsizes", samples)

    emit:
    samples
}