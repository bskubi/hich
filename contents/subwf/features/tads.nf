include {createCompositeStrategy; filterSamplesByStrategy} from '../util/analysisPlans.nf'
include {skip} from '../util/cli.nf'
include {withLog; stubLog} from '../util/logs.nf'

process HiCExplorerFindTADs {
    publishDir "results/tads",
               mode: params.general.publish.mode
    label 'features'
    container params.general.hicexplorerContainer
    tag "$id"

    input:
    tuple val(id), path(mcool), val(resolution), val(hicExplorerFindTADsParams)

    output:
    tuple val(id), path("${id}_boundaries.bed"), path("${id}_boundaries.gff"), path("${id}_domains.bed"), path("${id}_score.bedgraph"), path("${id}_tad_score.bm")

    shell:
    cmd = ["hicFindTADs --outPrefix ${id} --matrix ${mcool}::/resolutions/${resolution} -p ${task.cpus}"] + hicExplorerFindTADsParams
    cmd = cmd.join(" ")

    logMap = [
        task: "HiCExplorerFindTADs",
        input: [id: id, mcool: mcool, resolution: resolution, hicExplorerFindTADsParams: hicExplorerFindTADsParams], 
        output: [
            boundariesBed: "${id}_boundaries.bed",
            boundariesGff: "${id}_boundaries.gff",
            domainsBed: "${id}_domains.bed",
            scoreBedgraph: "${id}_score.bedgraph",
            tadScoreBm: "${id}_tad_score.bm"
        ]
    ]
    
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}_boundaries.bed' '${id}_boundaries.gff' '${id}_domains.bed' '${id}_score.bedgraph' '${id}_tad_score.bm'"
    cmd = ["hicFindTADs --outPrefix ${id} --matrix ${mcool}::/resolutions/${resolution} -p ${task.cpus}"] + hicExplorerFindTADsParams
    cmd = cmd.join(" ")

    logMap = [
        task: "HiCExplorerFindTADs",
        input: [id: id, mcool: mcool, resolution: resolution, hicExplorerFindTADsParams: hicExplorerFindTADsParams], 
        output: [
            boundariesBed: "${id}_boundaries.bed",
            boundariesGff: "${id}_boundaries.gff",
            domainsBed: "${id}_domains.bed",
            scoreBedgraph: "${id}_score.bedgraph",
            tadScoreBm: "${id}_tad_score.bm"
        ]
    ]
    stubLog(stub, cmd, logMap)
}

workflow TADs {
    take:
    samples

    main:

    if (!skip("tads")) {
        params.tads.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.mcool, analysisPlan.resolution, analysisPlan.hicExplorerFindTADsParams)
                }
                | HiCExplorerFindTADs
        }
    }


    emit:
    samples
}