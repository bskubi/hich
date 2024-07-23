include {MakeResourceFile} from "./makeResourceFile.nf"

process MakeDigest {
    //container "redigest"
    
    input:
    tuple path(reference), val(enzymes), val(fragfile), val(assembly)

    output:
    tuple path(reference), val(enzymes), path(fragfile), val(assembly)

    script:
    "redigest --output ${fragfile} ${reference} ${enzymes}"
}

workflow MakeMissingDigest {
    take:
        samples
    
    main:
        def digestEnzymesDeclared = {it.get("enzymes").trim().length() != 0}

        def hasFragfileName = {it.get("fragfile").trim().length() > 0}
        
        def fragfileExists = {hasFragfileName(it) && file(it.fragfile).exists()}

        def ensureFragfileName = {
            if (digestEnzymesDeclared(it) && !hasFragfileName(it)) {
                it.fragfile = "${it.assembly}_${it.enzymes}.bed"
            }
            it;
        }

        samples
            | map{ensureFragfileName(it)}
            | branch{
                exists: digestEnzymesDeclared(it) && fragfileExists(it)
                missing: digestEnzymesDeclared(it) && !fragfileExists(it)
                no_change: true
            }
            | set{branched}
        
        samples = MakeResourceFile(
            branched.exists,
            branched.missing,
            branched.no_change,
            "fragfile",
            MakeDigest,
            ["reference", "enzymes", "fragfile", "assembly"],
            ["reference", "enzymes", "fragfile", "assembly"],
            ["assembly", "enzymes"]
        )

    emit:
        samples
}