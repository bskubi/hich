include {QCReads} from './qcHicReads.nf'
include {transpack; emptyOnLastStep; updateChannel; pack2} from './extraops.nf'

process HicToMcool {
    publishDir params.general.publish.mcool ? params.general.publish.mcool : "results",
               saveAs: {params.general.publish.mcool ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::hictk"
    container "ghcr.io/paulsengroup/hictk:1.0.0"
    label 'doJobArray'
    label 'convertHicToMcool'

    input:
    tuple val(id), path(hicFile)

    output:
    tuple val(id), path(mcoolFile)

    shell:

    hicFileName = hicFile.getFileName().toString()
    prefix = hicFileName.substring(0, hicFileName.lastIndexOf("."))
    mcoolFile = "${prefix}.mcool"
    cpus = task.cpus >= 2 ? task.cpus : 2
    threads = task.cpus ? "-t ${cpus}" : ""
    "hictk convert ${threads} ${hicFile} ${mcoolFile}"

    stub:
    hicFileName = hicFile.getFileName().toString()
    prefix = hicFileName.substring(0, hicFileName.lastIndexOf("."))
    mcoolFile = "${prefix}.mcool"
    "touch ${hicFile} ${mcoolFile}"

}

process McoolToHic {
    publishDir params.general.publish.hic ? params.general.publish.hic : "results",
               saveAs: {params.general.publish.hic ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::hictk"
    container "ghcr.io/paulsengroup/hictk:1.0.0"
    label 'doJobArray'
    label 'convertMcoolToHic'

    input:
    tuple val(id), path(mcoolFile)

    output:
    tuple val(id), path(hicFile)

    shell:

    mcoolFileName = mcoolFile.getFileName().toString()
    prefix = mcoolFileName.substring(0, mcoolFileName.lastIndexOf("."))
    hicFile = "${prefix}.hic"
    print(task.cpus)
    cpus = task.cpus >= 2 ? task.cpus : 2
    threads = task.cpus ? "-t ${cpus}" : ""
    "hictk convert ${threads} ${mcoolFile} ${hicFile}"

    stub:
    mcoolFileName = mcoolFile.getFileName().toString()
    prefix = mcoolFileName.substring(0, mcoolFileName.lastIndexOf("."))
    hicFile = "${prefix}.hic"
    "touch ${hicFile}"

}

workflow IngestMatrix {
    take:
        samples

    main:

    samples
        | filter{it.datatype == "matrix" && it.mcool && !it.hic}
        | map{tuple(it.id, it.mcool)}
        | McoolToHic
        | map{
            id, hicFile ->
            [id:id, hic:hicFile, latest:hicFile, latestMatrix:hicFile]}
        | set{mcoolToHic}
    
    samples
        | filter{it.datatype == "matrix" && it.hic && !it.mcool}
        | map{tuple(it.id, it.hic)}
        | HicToMcool
        | map{
            id, mcoolFile ->
            [id:id, mcool:mcoolFile, latest:mcoolFile, latestMatrix:mcoolFile]}
        | set{hicToMcool}
    
    pack2(samples, mcoolToHic) | set{samples}
    pack2(samples, hicToMcool) | set{samples}

    samples = emptyOnLastStep("IngestMatrix", samples)

    emit:
        samples
}