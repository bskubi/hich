include {MakeResourceFile} from "./makeResourceFile.nf"
include {source} from "./extraops.nf"

process MakeDigest {
    container "redigest"
    
    input:
    tuple path(reference), val(enzymes), val(fragfile), val(assembly)

    output:
    tuple val(enzymes), path(fragfile), val(assembly)

    script:
    "redigest --output ${fragfile} ${reference} ${enzymes}"

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
           {"${it.assembly}_${it.enzymes}.bed"},
           [["assembly", "enzymes"]],
           {it.enzymes}) | set{samples}
    
    emit:
    samples
}