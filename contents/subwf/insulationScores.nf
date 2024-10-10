include {createCompositeStrategy; filterSamplesByStrategy; skip} from './extraops.nf'

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

workflow InsulationScores {
    take:
    samples

    main:

    if (!skip("insulationScores")) {
        params.insulation.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)
            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.latestMatrix, analysisPlan.resolution, analysisPlan.cooltoolsInsulationParams)
                }
                | CooltoolsInsulation
        }
    }


    emit:
    samples
}