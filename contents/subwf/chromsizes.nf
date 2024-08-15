include {source; emptyOnLastStep} from './extraops.nf'

process MakeChromsizes {
    publishDir params.general.publish.chromsizes ? params.general.publish.chromsizes : "results",
               saveAs: {params.general.publish.chromsizes ? it : null},
               mode: params.general.publish.mode

    conda "bioconda::ucsc-fasize"
    container 'quay.io/biocontainers/ucsc-fasize:332--0'

    input:
    tuple path(reference), val(assembly), val(chromsizes)

    output:
    tuple val(assembly), path(chromsizes)

    shell:
    "faSize -detailed -tab ${reference} > ${chromsizes}"

    stub:
    "touch ${chromsizes}"
}

workflow MakeMissingChromsizes {
    take:
    samples

    main:

    source(MakeChromsizes,
           samples,
           "chromsizes",
           ["genomeReference", "assembly", "chromsizes"],
           ["assembly", "chromsizes"],
           {"${it.assembly}.sizes"},
           "assembly",
            {true}) | set{samples}
    
    samples = emptyOnLastStep("chromsizes", samples)

    emit:
    samples
}