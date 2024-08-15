include {source; emptyOnLastStep} from "./extraops.nf"

process FragmentIndexProc {
    publishDir params.general.publish.fragmentIndex ? params.general.publish.fragmentIndex : "results",
               saveAs: {params.general.publish.fragmentIndex ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"
    
    input:
    tuple path(genomeReference), val(restriction_enzymes), val(fragfile), val(assembly)

    output:
    tuple val(restriction_enzymes), path(fragfile), val(assembly)

    script:
    "hich digest --output ${fragfile} ${genomeReference} ${restriction_enzymes}"

    stub:
    "touch ${fragfile}"
}

workflow FragmentIndex {
    take:
    samples
    
    main:
    
    source(FragmentIndexProc,
           samples,
           "fragfile",
            ["genomeReference", "restriction_enzymes", "fragfile", "assembly"],
            ["restriction_enzymes", "fragfile", "assembly"],
           {"${it.assembly}_${it.restriction_enzymes.replace(" ", "_")}.bed"},
           [["assembly", "restriction_enzymes"]],
           {it.restriction_enzymes}) | set{samples}
    
    samples = emptyOnLastStep("fragmentIndex", samples)
    
    emit:
    samples
}