include {emptyOnLastStep; pack; skip} from './extraops.nf'

process JuicerToolsPre {
    /*
        Juicer Tools Pre documentation: https://github.com/aidenlab/juicer/wiki/Pre

        We use version 1.22.01 as 2.0+ versions are in development and certain
        features available in version 1 are unavailable in 2.
    */
    publishDir params.general.publish.hic ? params.general.publish.hic : "results",
               saveAs: {params.general.publish.hic ? it : null},
               mode: params.general.publish.mode
    container params.general.juicerContainer
    label 'createMatrix'
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
            bins] + juicerToolsPreParams + ["'${infile}' '${outfile}' '${chromsizes}'"]
    cmd.removeAll([null])
    cmd.join(" ")

    stub:
    "touch ${id}.hic"
}

workflow HicMatrix {
    take:
    samples
    
    main:
    samples
        | filter{!skip("hicMatrix") && it.matrix.makeHicFileFormat && (it.pairs || it.latestPairs) && !it.hic}
        | map{tuple(it.id, it.latest, it.chromsizes, it.pairsFormat, it.matrix, it.juicerToolsPreParams)}
        | JuicerToolsPre
        | map{id, hic -> [id: id, hic: hic, latestMatrix: hic]}
        | set{result}
    pack(samples, result) | set{samples}

    samples = emptyOnLastStep("hicMatrix", samples)

    emit:
    samples
}