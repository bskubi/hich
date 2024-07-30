process CooltoolsInsulation {
    publishDir "results/insulation"
    container "bskubi/open2c:latest"

    input:
    tuple val(id), path(mcool), val(resolution), val(cooltools_insulation_params)

    output:
    tuple val(id), path("${id}_insulation.tsv"), path("${id}_insulation.tsv.${resolution}.bw")

    shell:
    cmd = ["cooltools insulation", 
           "--output ${id}_insulation.tsv"] +
          cooltools_insulation_params +
          ["${mcool}::/resolutions/${resolution} 100000"]
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${id}_insulation.tsv ${id}_insulation.tsv.${resolution}.bw"
}

workflow CallInsulation {
    take:
    samples

    main:
    samples
        | filter {it.get("insulation") != null && it.get("mcool") != null}
        | map{tuple(it.id, it.mcool, it.insulation.resolution, it.insulation.cooltools_insulation_params)}
        | CooltoolsInsulation

    emit:
    samples
}