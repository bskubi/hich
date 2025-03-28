include {withLog; stubLog; emptyOnLastStep; pack; isExistingFile; skip} from "../extraops.nf"

process HichDigest {
    publishDir params.general.publish.fragmentIndex ? params.general.publish.fragmentIndex : "results",
               saveAs: {params.general.publish.fragmentIndex ? it : null},
               mode: params.general.publish.mode
    container params.general.hichContainer
    label 'smallResource'
    
    input:
    tuple val(genomeReferenceString), path(genomeReference), val(restrictionEnzymes), val(fragmentIndex), val(assembly)

    output:
    tuple val(genomeReferenceString), val(restrictionEnzymes), path(fragmentIndex), val(assembly)

    script:
    cmd = "hich digest --output '${fragmentIndex}' '${genomeReference}' ${restrictionEnzymes.split(',').join(' ')}"
    logMap = [task: "HichDigest", input: [genomeReference: genomeReference, restrictionEnzymes: restrictionEnzymes, fragmentIndex: fragmentIndex, assembly: assembly], 
    output: [fragmentIndex: fragmentIndex]]
    withLog(cmd, logMap)

    stub:
    stub = "touch '${fragmentIndex}'"
    cmd = "hich digest --output '${fragmentIndex}' '${genomeReference}' ${restrictionEnzymes.split(',').join(' ')}"
    logMap = [task: "HichDigest", input: [genomeReference: genomeReference, restrictionEnzymes: restrictionEnzymes, fragmentIndex: fragmentIndex, assembly: assembly], 
    output: [fragmentIndex: fragmentIndex]]
    stubLog(stub, cmd, logMap)
}

workflow FragmentIndex {
    take:
    samples
    
    main:
    
    if (!skip("fragmentIndex")) {
        samples
            | filter{it.restrictionEnzymes && !isExistingFile(it.fragmentIndex)}
            | map{it.fragmentIndex = it.fragmentIndex ?: "${it.assembly}_${it.restrictionEnzymes.replace(" ", "_")}.bed"; it}
            | map{tuple(it.genomeReference, it.genomeReference, it.restrictionEnzymes, it.fragmentIndex, it.assembly)}
            | unique
            | HichDigest
            | map{genomeReference, restrictionEnzymes, fragmentIndex, assembly -> 
                [genomeReference: file(genomeReference), restrictionEnzymes: restrictionEnzymes, fragmentIndex: fragmentIndex, assembly: assembly]}
            | set{result}
        pack(samples, result, ["genomeReference", "assembly", "restrictionEnzymes"]) | set{samples}
    }

    
    samples = emptyOnLastStep("fragmentIndex", samples)
    
    emit:
    samples
}