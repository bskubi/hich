include {source} from "./extraops.nf"

process MakeDigest {
    publishDir params.general.publish.digest ? params.general.publish.digest : "results",
               saveAs: {params.general.publish.digest ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"
    
    input:
    tuple path(reference), val(enzymes), val(fragfile), val(assembly)

    output:
    tuple val(enzymes), path(fragfile), val(assembly)

    script:
    "hich digest --output ${fragfile} ${reference} ${enzymes}"

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
            ["reference", "enzymes", "fragfile", "assembly"],
            ["enzymes", "fragfile", "assembly"],
           {"${it.assembly}_${it.enzymes.replace(" ", "_")}.bed"},
           [["assembly", "enzymes"]],
           {it.enzymes}) | set{samples}
    
    emit:
    samples
}