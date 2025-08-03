include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {withLog; stubLog} from '../util/logs.nf'

process HicToMcool {
    publishDir params.general.publish.mcool ? params.general.publish.mcool : "results",
               saveAs: {params.general.publish.mcool ? it : null},
               mode: params.general.publish.mode
    
    label 'convertHicToMcool'
    cpus 2
    maxForks 1  // Necessary due to unexplained bug in hictk
    conda "$projectDir/env/dev_env.yml"
    container params.general.hictkContainer

    input:
    tuple val(id), path(hicFile)

    output:
    tuple val(id), path(mcoolFile)

    shell:

    hicFileName = hicFile.getFileName().toString()
    prefix = hicFileName.substring(0, hicFileName.lastIndexOf("."))
    mcoolFile = "${prefix}.mcool"
    cmd = "mkdir ./tmp && hictk convert --tmpfile ./tmp '${hicFile}' '${mcoolFile}'"
    logMap = [task: "HicToMcool", input: [id: id, hicFile: hicFile], 
    output: [mcoolFile: mcoolFile]]
    withLog(cmd, logMap)

    stub:
    hicFileName = hicFile.getFileName().toString()
    prefix = hicFileName.substring(0, hicFileName.lastIndexOf("."))
    mcoolFile = "${prefix}.mcool"
    stub = "touch '${hicFile}' '${mcoolFile}'"
    cmd = "mkdir ./tmp && hictk convert --tmpfile ./tmp '${hicFile}' '${mcoolFile}'"
    logMap = [task: "HicToMcool", input: [id: id, hicFile: hicFile], 
    output: [mcoolFile: mcoolFile]]
    stubLog(stub, cmd, logMap)

}

process McoolToHic {
    publishDir params.general.publish.hic ? params.general.publish.hic : "results",
               saveAs: {params.general.publish.hic ? it : null},
               mode: params.general.publish.mode
    
    label 'convertMcoolToHic'
    conda "$projectDir/env/dev_env.yml"
    container params.general.hictkContainer

    input:
    tuple val(id), path(mcoolFile)

    output:
    tuple val(id), path(hicFile)

    shell:

    mcoolFileName = mcoolFile.getFileName().toString()
    prefix = mcoolFileName.substring(0, mcoolFileName.lastIndexOf("."))
    hicFile = "${prefix}.hic"
    cpus = task.cpus && task.cpus >= 2 ? task.cpus : 2
    threads = cpus > 2 ? "-t ${cpus}" : ""
    cmd = "mkdir ./tmp && hictk convert --tmpfile ./tmp ${threads} '${mcoolFile}' '${hicFile}'"
    logMap = [task: "McoolToHic", input: [id: id, mcoolFile: mcoolFile], 
    output: [hicFile: hicFile]]
    withLog(cmd, logMap)

    stub:
    mcoolFileName = mcoolFile.getFileName().toString()
    prefix = mcoolFileName.substring(0, mcoolFileName.lastIndexOf("."))
    hicFile = "${prefix}.hic"
    stub = "touch '${hicFile}'"
    cpus = task.cpus && task.cpus >= 2 ? task.cpus : 2
    threads = cpus > 2 ? "-t ${cpus}" : ""
    cmd = "mkdir ./tmp && hictk convert --tmpfile ./tmp ${threads} '${mcoolFile}' '${hicFile}'"
    logMap = [task: "McoolToHic", input: [id: id, mcoolFile: mcoolFile], 
    output: [hicFile: hicFile]]
    stubLog(stub, cmd, logMap)

}

workflow IngestMatrix {
    take:
        samples

    main:

    if (!skip("ingestMatrix")) {
        samples
            | filter{it.datatype == "matrix" && makeHicFileFormat && it.mcool && !it.hic}
            | map{tuple(it.id, it.mcool)}
            | McoolToHic
            | map{
                id, hicFile ->
                [id:id, hic:hicFile, latest:hicFile, latestMatrix:hicFile]}
            | set{mcoolToHic}
        
        samples
            | filter{it.datatype == "matrix" && makeMcoolFileFormat && it.hic && !it.mcool}
            | map{tuple(it.id, it.hic)}
            | HicToMcool
            | map{
                id, mcoolFile ->
                [id:id, mcool:mcoolFile, latest:mcoolFile, latestMatrix:mcoolFile]}
            | set{hicToMcool}

        keyUpdate(samples, mcoolToHic, "id") | set{samples}
        keyUpdate(samples, hicToMcool, "id") | set{samples}
    }

    


    samples = emptyOnLastStep("ingestMatrix", samples)

    emit:
        samples
}