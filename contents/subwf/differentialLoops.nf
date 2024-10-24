include {createCompositeStrategy; pairSamplesByStrategy; skip} from './extraops.nf'

process MustacheDiffloops{
    publishDir "results/loops",
               mode: params.general.publish.mode
    container params.general.mustacheContainer
    label 'features'

    input:
    tuple val(prefix), path(mx1, stageAs: 'mx1/*'), path(mx2, stageAs: 'mx2/*'), val(mustacheParams)

    output:
    tuple path("${prefix}.loop1"), path("${prefix}.loop2"), path("${prefix}.diffloop1"), path("${prefix}.diffloop2")

    shell:
    cmd = ["python /mustache/mustache/diff_mustache.py",
           "-f1 ${mx1} -f2 ${mx2}",
           "-o ${prefix}"] + mustacheParams
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch ${prefix}.loop1 ${prefix}.loop2 ${prefix}.diffloop1 ${prefix}.diffloop2"
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
                | map{
                    s1, s2 ->
                    prefix = "${s1.id}_${s2.id}"
                    tuple(prefix, s1.latestMatrix, s2.latestMatrix, analysisPlan.mustacheParams)
                }
                | MustacheDiffloops
        }
    }


    emit:
    samples
}