process HichCompartments {
    publishDir "results/compartments",
               mode: params.general.publish.mode

    container "bskubi/hich:latest"

    input:
    tuple val(id), path(reference), path(matrix), val(resolution), val(hich_compartments_params)

    output:
    tuple val(id), path("${id}_0.bw"), path("${id}_1.bw"), path("${id}_2.bw")

    shell:
    cmd = ["hich compartments"] +
          hich_compartments_params +
          ["${reference} ${matrix} ${resolution}"]
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${id}_0.bw ${id}_1.bw ${id}_2.bw"
}

workflow CallCompartments {
    take:
    samples

    main:
    samples
        | filter {it.get("compartments") != null && it.get("latest_matrix") != null}
        | map{tuple(it.id, it.reference, it.latest_matrix, it.compartments.resolution, it.compartments.hich_compartments_params)}
        | HichCompartments

    emit:
    samples
}
