include {MakeResourceFile} from './makeResourceFile.nf'

process MakeChromsizes {
    conda 'bioconda::ucsc-fasize'

    input:
    tuple path(reference), val(assembly), val(chromsizes)

    output:
    tuple val(assembly), path(chromsizes)

    shell:
    "faSize -detailed -tab ${reference} > ${chromsizes}"
}

workflow MakeMissingChromsizes {
    take:
        samples
    
    main:
        def hasChromsizesFilename = {
            it.get("chromsizes", "").toString().trim().length() > 0
            && it.get("chromsizes") != "NULL"}

        def chromsizesExists = {hasChromsizesFilename(it)
                                && file(it.chromsizes).exists()}
        
        def ensureChromsizesFilename = {
            if (!hasChromsizesFilename(it)) {
                it.chromsizes = "${it.assembly}.sizes"
            }
            it
        }

        samples
            | map{ensureChromsizesFilename(it)}
            | branch{
                exists: chromsizesExists(it)
                missing: !chromsizesExists(it) && hasChromsizesFilename(it)
                no_change: true
            }
            | set{branched}
        
        samples = MakeResourceFile(
            branched.exists,
            branched.missing,
            branched.no_change,
            "chromsizes",
            MakeChromsizes,
            ["reference", "assembly", "chromsizes"],
            ["assembly", "chromsizes"],
            ["assembly"]
        )

    emit:
        samples
}