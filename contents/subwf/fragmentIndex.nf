include {source; emptyOnLastStep; pack2; isExistingFile} from "./extraops.nf"

process FragmentIndexProc {
    publishDir params.general.publish.fragmentIndex ? params.general.publish.fragmentIndex : "results",
               saveAs: {params.general.publish.fragmentIndex ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"
    label 'smallResource'
    
    input:
    tuple path(genomeReference), val(restrictionEnzymes), val(fragmentIndex), val(assembly)

    output:
    tuple val(restrictionEnzymes), path(fragmentIndex), val(assembly)

    script:
    "hich digest --output ${fragmentIndex} ${genomeReference} ${restrictionEnzymes}"

    stub:
    "touch ${fragmentIndex}"
}

workflow FragmentIndex {
    take:
    samples
    
    main:
    
    samples
        | filter{it.restrictionEnzymes && !it.isExistingFile(it.fragmentIndex)}
        | map{it.fragmentIndex = it.fragmentIndex ?: "${it.assembly}_${it.restrictionEnzymes.replace(" ", "_")}.bed"; it}
        | map{tuple(it.genomeReference, it.restrictionEnzymes, it.fragmentIndex, it.assembly)}
        | unique
        | FragmentIndexProc
        | map{restrictionEnzymes, fragmentIndex, assembly -> 
              [restrictionEnzymes: restrictionEnzymes, fragmentIndex: fragmentIndex, assembly: assembly]}
        | set{result}
    pack2(samples, result, ["assembly", "restrictionEnzymes"]) | set{samples}
    
    samples = emptyOnLastStep("FragmentIndex", samples)
    
    emit:
    samples
}