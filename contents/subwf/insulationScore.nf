process CooltoolsInsulation {
    publishDir "results/insulation",
               mode: params.general.publish.mode
    container "bskubi/hich:latest"
    conda "bioconda::cooltools"

    input:
    tuple val(id), path(mcool), val(resolution), val(cooltoolsInsulationParams)

    output:
    tuple val(id), path("${id}_insulation.tsv"), path("${id}_insulation.tsv.${resolution}.bw")

    shell:
    cmd = ["cooltools insulation", 
           "--output ${id}_insulation.tsv"] +
          cooltoolsInsulationParams +
          ["${mcool}::/resolutions/${resolution} 100000"]
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${id}_insulation.tsv ${id}_insulation.tsv.${resolution}.bw"
}

workflow InsulationScore {
    take:
    samples

    main:
    samples
        | filter {it.get("insulation") != null && it.get("mcool") != null}
        | map{tuple(it.id, it.mcool, it.insulation.resolution, it.insulation.cooltoolsInsulationParams)}
        | CooltoolsInsulation

    emit:
    samples
}