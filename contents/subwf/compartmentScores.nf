include {createCompositeStrategy; filterSamplesByStrategy; skip} from './extraops.nf'


process HichCompartments {
    publishDir "results/compartments",
               mode: params.general.publish.mode

    container params.general.hichContainer
    label 'features'


    input:
    tuple val(id), path(genomeReference), path(matrix), val(resolution), val(hichCompartmentsParams)

    output:
    tuple val(id), path("${id}_0.bw"), path("${id}_1.bw"), path("${id}_2.bw")

    shell:
    cmd = ["hich compartments --n_eigs 3"] +
          hichCompartmentsParams +
          ["${genomeReference} ${matrix} ${resolution}"]
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${id}_0.bw ${id}_1.bw ${id}_2.bw"
}

workflow CompartmentScores {
    take:
    samples

    main:
    if (!skip("compartmentScores")) {
        params.compartments.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.genomeReference, sample.latestMatrix, analysisPlan.resolution, analysisPlan.hichCompartmentsParams)
                }
                | HichCompartments
        }
    }

    emit:
    samples
}
