include {createCompositeStrategy; filterSamplesByStrategy;} from '../util/analysisPlans.nf'
include {skip} from '../util/cli.nf'
include {withLog; stubLog} from '../util/logs.nf'

process MustacheLoops{
    publishDir "results/loops",
               mode: params.general.publish.mode
    container params.general.mustacheContainer
    label 'features'
    tag "$id"

    input:
    tuple val(id), val(prefix), path(mx), val(mustacheParams)

    output:
    path("${prefix}.loop")

    shell:
    cmd = ["mustache -f '${mx}' -o '${prefix}.loop'"] + mustacheParams
    cmd = cmd.join(" ")
    logMap = [task: "MustacheLoops", input: [id: id, prefix: prefix, mx: mx, mustacheParams: mustacheParams], 
    output: [loops: "${prefix}.loop"]]
    withLog(cmd, logMap)

    stub:
    stub = "touch '${prefix}.loop'"
    cmd = ["mustache -f '${mx}' -o '${prefix}.loop'"] + mustacheParams
    cmd = cmd.join(" ")
    logMap = [task: "MustacheLoops", input: [id: id, prefix: prefix, mx: mx, mustacheParams: mustacheParams], 
    output: [loops: "${prefix}.loop"]]
    stubLog(stub, cmd, logMap)
}

workflow Loops {
    take:
    samples

    main:
    
    if (!skip("loops")) {
        params.loops.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.id, sample.latestMatrix, analysisPlan.mustacheParams)
                }
                | MustacheLoops
        }
    }


    emit:
    samples
}