include {createCompositeStrategy; pairSamplesByStrategy} from '../util/analysisPlans.nf'
include {skip} from '../util/cli.nf'
include {withLog; stubLog} from '../util/logs.nf'

process MustacheDiffloops{
    publishDir "results/loops",
               mode: params.general.publish.mode
    container params.general.mustacheContainer
    label 'features'
    tag "$id1 $id2"

    input:
    tuple val(id1), val(id2), val(prefix), path(mx1, stageAs: 'mx1/*'), path(mx2, stageAs: 'mx2/*'), val(mustacheParams)

    output:
    tuple path("${prefix}.loop1"), path("${prefix}.loop2"), path("${prefix}.diffloop1"), path("${prefix}.diffloop2")

    shell:
    cmd = ["diff_mustache",
           "-f1 '${mx1}' -f2 '${mx2}'",
           "-o '${prefix}'"] + mustacheParams
    cmd = cmd.join(" ")
    logMap = [task: "MustacheDiffloops", input: [id1: id1, id2: id2, mx1: mx1, mx2: mx2, mustacheParams: mustacheParams], output: [loop1: "${prefix}.loop1", loop2: "${prefix}.loop2", diffloop1: "${prefix}.diffloop1", diffloop2: "${prefix}.diffloop2"]]
    withLog(cmd, logMap)

    stub:
    stub = "touch '${prefix}.loop1' '${prefix}.loop2' '${prefix}.diffloop1' '${prefix}.diffloop2'"
    cmd = ["diff_mustache",
           "-f1 '${mx1}' -f2 '${mx2}'",
           "-o '${prefix}'"] + mustacheParams
    cmd = cmd.join(" ")
    logMap = [task: "MustacheDiffloops", input: [id1: id1, id2: id2, mx1: mx1, mx2: mx2, mustacheParams: mustacheParams], output: [loop1: "${prefix}.loop1", loop2: "${prefix}.loop2", diffloop1: "${prefix}.diffloop1", diffloop2: "${prefix}.diffloop2"]]
    stubLog(stub, cmd, logMap)

}

workflow DifferentialLoops {
    take:
    samples

    main:
    
    if (!skip("differentialLoops")) {
        params.differentialLoops.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            pairSamplesByStrategy(samples, strategy)
                | filter{s1, s2 -> s1.id != s2.id && s1.matrixPlanName == s2.matrixPlanName}
                | map{
                    s1, s2 ->
                    prefix = "${s1.id}_${s2.id}"
                    tuple(s1.id, s2.id, prefix, s1.latestMatrix, s2.latestMatrix, analysisPlan.mustacheParams)
                }
                | MustacheDiffloops
        }
    }


    emit:
    samples
}