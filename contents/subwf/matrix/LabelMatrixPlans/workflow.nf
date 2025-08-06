include {createCompositeStrategy; filterSamplesByStrategy} from '../../util/analysisPlans.nf'
include {columnsToRows} from '../../util/reshape.nf'
include {emptyOnLastStep} from '../../util/cli.nf'
include {Setup} from '../../setup/setup.nf'

workflow LabelMatrixPlans {
    take:
    samples

    main:
    newSamples = channel.empty()

    params.matrices.each {
        planName, analysisPlan ->
        analysisPlan["plan_name"] = planName
        strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

        filterSamplesByStrategy(samples, strategy)
            | map{it + [matrix_opts: analysisPlan, id: it.id + "_${planName}"]}
            | concat(newSamples)
            | set{newSamples}
        
    }
    // TODO: Warn about no matrices.
    if (!params.containsKey("keepNoMatrixPlan")) {
        newSamples
            | set{samples}
    } else {
        newSamples
            | concat(samples)
            | set{samples}
    }

    emit:
    samples
}