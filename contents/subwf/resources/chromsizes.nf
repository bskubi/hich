include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {isExistingFile} from '../util/files.nf'
include {withLog; stubLog} from '../util/logs.nf'

process FASIZE {
    publishDir params.general.publish.chromsizes ? params.general.publish.chromsizes : "results",
               saveAs: {params.general.publish.chromsizes ? it : null},
               mode: params.general.publish.mode

    container params.general.chromsizesContainer
    label 'smallResource'

    input:
    tuple val(genomeReferenceString), path(genomeReference), val(assembly), val(chromsizes)

    output:
    tuple val(genomeReferenceString), val(assembly), path(chromsizes)

    shell:
    cmd = "faSize -detailed -tab '${genomeReference}' > '${chromsizes}'"
    logMap = [task: "FaSize", "input": [genomeReference: genomeReference, assembly: assembly], "output": chromsizes].collectEntries { k, v -> [k, v.toString()] }
    withLog(cmd, logMap)

    stub:
    stub = "touch '${chromsizes}'"
    cmd = "faSize -detailed -tab '${genomeReference}' > '${chromsizes}'"

    logMap = [task: "FaSize", "input": [genomeReference: genomeReference, assembly: assembly], "output": chromsizes].collectEntries { k, v -> [k, v.toString()] }
    stubLog(stub, cmd, logMap)
}

workflow Chromsizes {
    take:
    samples

    main:
    
    if (!skip("chromsizes")) {
        samples
            | filter {!isExistingFile(it.chromsizes)}
            | map{tuple(it.genomeReference, it.genomeReference, it.assembly, "${it.assembly}.sizes")}
            | unique
            | FASIZE
            | map{genomeReference, assembly, chromsizes -> [genomeReference: file(genomeReference), assembly: assembly, chromsizes: chromsizes]}
            | set{result}
        keyUpdate(samples, result, ["genomeReference", "assembly"]) | set{samples}
    }


    samples = emptyOnLastStep("chromsizes", samples)

    emit:
    samples
}