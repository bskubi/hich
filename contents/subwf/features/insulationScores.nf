include {createCompositeStrategy; filterSamplesByStrategy; skip} from '../extraops.nf'
include {withLog; stubLog} from '../util/logs.nf'

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
    cmd = ["cooltools insulation --bigwig", 
           "--output '${id}_insulation.tsv'"] +
          cooltoolsInsulationParams +
          ["'${mcool}::/resolutions/${resolution}' ${window}"]
    cmd = cmd.join(" ")
    logMap = [task: "CooltoolsInsulation", input: [id: id, mcool: mcool, resolution: resolution, cooltoolsInsulationParams: cooltoolsInsulationParams, window: window], 
    output: [insulationTSV: "${id}_insulation.tsv", insulationBW: "${id}_insulation.tsv.${resolution}.bw"]]
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}_insulation.tsv' '${id}_insulation.tsv.${resolution}.bw'"
    cmd = ["cooltools insulation", 
           "--output '${id}_insulation.tsv'"] +
          cooltoolsInsulationParams +
          ["'${mcool}::/resolutions/${resolution}' ${window}"]
    cmd = cmd.join(" ")
    logMap = [task: "CooltoolsInsulation", input: [id: id, mcool: mcool, resolution: resolution, cooltoolsInsulationParams: cooltoolsInsulationParams, window: window], 
    output: [insulationTSV: "${id}_insulation.tsv", insulationBW: "${id}_insulation.tsv.${resolution}.bw"]]
    stubLog(stub, cmd, logMap)
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
                    tuple(sample.id, sample.latestMatrix, analysisPlan.resolution, analysisPlan.cooltoolsInsulationParams, analysisPlan.window)
                }
                | CooltoolsInsulation
        }
    }


    emit:
    samples
}