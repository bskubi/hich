include {transpack; emptyOnLastStep} from './extraops.nf'

process JuicerToolsPre {
    /*
        Juicer Tools Pre documentation: https://github.com/aidenlab/juicer/wiki/Pre

        We use version 1.22.01 as 2.0+ versions are in development and certain
        features available in version 1 are unavailable in 2.
    */
    publishDir params.general.publish.hic ? params.general.publish.hic : "results",
               saveAs: {params.general.publish.hic ? it : null},
               mode: params.general.publish.mode
    container "bskubi/juicer_tools:1.22.01"
    cpus 5
    memory {2.GB * task.attempt}

    input:
    tuple val(id), path(infile), path(chromsizes), val(pairsFormat), val(matrix), val(juicerToolsPreParams)

    output:
    tuple val(id), path("${id}.hic")

    shell:
    outfile = "${id}.hic"
    min_bin = matrix.resolutions.min()
    bins = matrix.resolutions ? "-r ${matrix.resolutions.join(',')}" : ""

    cmd = ["java -Xmx2g -jar /app/juicer_tools_1.22.01.jar pre",
            bins] + juicerToolsPreParams + ["${infile} ${outfile} ${chromsizes}"]
    cmd.removeAll([null])
    cmd.join(" ")

    stub:
    "touch ${id}.hic"
}

workflow HicMatrix {
    take:
    samples
    
    main:
    samples | filter{it.matrix.makeHicFileFormat} | set{hic}

    transpack (
        JuicerToolsPre,
        [hic, samples],
        ["id", "latest", "chromsizes", "pairsFormat", "matrix", "juicerToolsPreParams"],
        ["id", "hic"],
        ["latestMatrix":"hic"],
        "id",
        ["nullOk":"juicerToolsPreParams"]
    ) | set{samples}

    samples = emptyOnLastStep("HicMatrix", samples)

    emit:
    samples
}