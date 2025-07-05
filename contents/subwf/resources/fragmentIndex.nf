include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {isExistingFile} from '../util/files.nf'
include {withLog; stubLog} from '../util/logs.nf'
include {asList; cmdList} from '../util/dataStructures.nf'

process HichDigest {
    publishDir params.general.publish.fragmentIndex ? params.general.publish.fragmentIndex : "results",
               saveAs: {params.general.publish.fragmentIndex ? it : null},
               mode: params.general.publish.mode
    label 'smallResource'

    input:
    tuple val(genomeReferenceString), path(genomeReference), val(restrictionEnzymes), val(fragmentIndex), val(assembly)

    output:
    tuple val(genomeReferenceString), val(restrictionEnzymes), path(fragmentIndex), val(assembly)

    script:
    //cmd = "hich digest --output '${fragmentIndex}' '${genomeReference}' ${restrictionEnzymes.split(',').join(' ')}"
    enzymes = cmdList(restrictionEnzymes.split(',').toList())

    cmd = "hich fasta re-fragments '${genomeReference}' '${fragmentIndex}' ${enzymes}"
    logMap = [task: "HichDigest", input: [genomeReference: genomeReference, restrictionEnzymes: restrictionEnzymes, fragmentIndex: fragmentIndex, assembly: assembly], 
    output: [fragmentIndex: fragmentIndex]]
    withLog(cmd, logMap)

    stub:
    stub = "touch '${fragmentIndex}'"
    enzymes = cmdList(restrictionEnzymes.split(','))

    cmd = "hich fasta re-fragments '${genomeReference}' '${fragmentIndex}' ${enzymes}"
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
        keyUpdate(samples, result, ["genomeReference", "assembly", "restrictionEnzymes"]) | set{samples}
    }

    
    samples = emptyOnLastStep("fragmentIndex", samples)
    
    emit:
    samples
}