include {createCompositeStrategy; filterSamplesByStrategy; skip} from './extraops.nf'

process CooltoolsInsulation {
    publishDir "results/insulation",
               mode: params.general.publish.mode
    container params.general.hichContainer
    conda "bioconda::cooltools"
    label 'features'

    input:
    tuple val(id), path(mcool), val(resolution), val(cooltoolsInsulationParams), val(window)

    output:
    tuple val(id), path("${id}_insulation.tsv"), path("${id}_insulation.tsv.${resolution}.bw")

    shell:
    cmd = ["cooltools insulation", 
           "--output ${id}_insulation.tsv"] +
          cooltoolsInsulationParams +
          ["${mcool}::/resolutions/${resolution} ${window}"]
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${id}_insulation.tsv ${id}_insulation.tsv.${resolution}.bw"
}

workflow InsulationScores {
    take:
    samples

    main:

    if (!skip("insulationScores")) {
        params.insulation.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            strategy += ["same": [], "different": []]

            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.latestMatrix, analysisPlan.resolution, analysisPlan.cooltoolsInsulationParams, analysisPlan.window)
                }
                | CooltoolsInsulation
        }
    }


    emit:
    samples
}