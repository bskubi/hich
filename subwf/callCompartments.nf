process CooltoolsEigsCis {
    publishDir "results/compartments"
    container "bskubi/open2c:latest"

    input:
    tuple val(id), path(mcool), val(resolution), val(cooltools_ciseigs_params)

    output:
    tuple val(id), path("${id}.cis.bw"), path("${id}.cis.lam.txt"), path("${id}.cis.vecs.tsv")

    shell:
    cmd = ["cooltools eigs-cis", 
           "--out-prefix ${id}"] +
          cooltools_ciseigs_params +
          ["${mcool}::/resolutions/${resolution}"]
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${id}.cis.bw ${id}.cis.lam.txt ${id}.cis.vecs.tsv"
}

workflow CallCompartments {
    take:
    samples

    main:
    samples
        | filter {it.get("compartments") != null && it.get("mcool") != null}
        | map{tuple(it.id, it.mcool, it.compartments.resolution, it.compartments.cooltools_eigs_cis_params)}
        | CooltoolsEigsCis

    emit:
    samples
}
