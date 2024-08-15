include {source; emptyOnLastStep} from "./extraops.nf"

process FragmentIndexProc {
    publishDir params.general.publish.fragmentIndex ? params.general.publish.fragmentIndex : "results",
               saveAs: {params.general.publish.fragmentIndex ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"
    
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
    
    source(FragmentIndexProc,
           samples,
           "fragmentIndex",
            ["genomeReference", "restrictionEnzymes", "fragmentIndex", "assembly"],
            ["restriction_enzymes", "fragmentIndex", "assembly"],
           {"${it.assembly}_${it.restrictionEnzymes.replace(" ", "_")}.bed"},
           [["assembly", "restrictionEnzymes"]],
           {it.restrictionEnzymes}) | set{samples}
    
    samples = emptyOnLastStep("fragmentIndex", samples)
    
    emit:
    samples
}