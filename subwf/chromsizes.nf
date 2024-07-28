include {MakeResourceFile} from './makeResourceFile.nf'
include {source} from './extraops.nf'

process MakeChromsizes {
    conda 'bioconda::ucsc-fasize'

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
           ["reference", "assembly", "chromsizes"],
           ["assembly", "chromsizes"],
           {"${it.assembly}.sizes"},
           "assembly") | set{result}

    emit:
    result
}