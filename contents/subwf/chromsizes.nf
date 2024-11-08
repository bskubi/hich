include {emptyOnLastStep; pack; skip; isExistingFile; withLog; stubLog} from './extraops.nf'

process FaSize {
    publishDir params.general.publish.chromsizes ? params.general.publish.chromsizes : "results",
               saveAs: {params.general.publish.chromsizes ? it : null},
               mode: params.general.publish.mode

    container params.general.chromsizesContainer
    label 'smallResource'
    memory 8.GB

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
    
    samples
        | filter {!skip("chromsizes") && !isExistingFile(it.chromsizes)}
        | map{tuple(it.genomeReference, it.genomeReference, it.assembly, "${it.assembly}.sizes")}
        | unique
        | FaSize
        | map{genomeReference, assembly, chromsizes -> [genomeReference: file(genomeReference), assembly: assembly, chromsizes: chromsizes]}
        | set{result}
    pack(samples, result, ["genomeReference", "assembly"]) | set{samples}

    samples = emptyOnLastStep("chromsizes", samples)

    emit:
    samples
}