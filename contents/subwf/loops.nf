include {createCompositeStrategy; filterSamplesByStrategy; skip} from './extraops.nf'

process MustacheLoops{
    publishDir "results/loops",
               mode: params.general.publish.mode
    container params.general.mustacheContainer
    label 'features'

    input:
    tuple val(prefix), path(mx), val(mustacheParams)

    output:
    path("${prefix}.loop")

    shell:
    cmd = ["python -m mustache -f '${mx}' -o '${prefix}.loop'"] + mustacheParams
    cmd = cmd.join(" ")
    cmd

    stub:
    "touch '${prefix}.loop'"
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
                    tuple(sample.id, sample.latestMatrix, analysisPlan.mustacheParams)
                }
                | MustacheLoops
        }
    }


    emit:
    samples
}