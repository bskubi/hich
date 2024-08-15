include {source; emptyOnLastStep} from "./extraops.nf"

process MakeDigest {
    publishDir params.general.publish.digest ? params.general.publish.digest : "results",
               saveAs: {params.general.publish.digest ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"
    
    input:
    tuple path(reference), val(restriction_enzymes), val(fragfile), val(assembly)

    output:
    tuple val(restriction_enzymes), path(fragfile), val(assembly)

    script:
    "hich digest --output ${fragfile} ${reference} ${restriction_enzymes}"

    stub:
    "touch ${fragfile}"
}

workflow MakeMissingDigest {
    take:
    samples
    
    main:
    
    source(MakeDigest,
           samples,
           "fragfile",
            ["reference", "restriction_enzymes", "fragfile", "assembly"],
            ["restriction_enzymes", "fragfile", "assembly"],
           {"${it.assembly}_${it.restriction_enzymes.replace(" ", "_")}.bed"},
           [["assembly", "restriction_enzymes"]],
           {it.restriction_enzymes}) | set{samples}
    
    samples = emptyOnLastStep("fragfile") ?: samples
    
    emit:
    samples
}