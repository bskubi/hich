include {createCompositeStrategy; filterSamplesByStrategy} from '../../util/analysisPlans.nf'
include {skip; emptyOnLastStep} from '../../util/cli.nf'
include {TADS} from './process.nf'

workflow TADs {
    take:
    samples

    main:

    myName = "TADs"
    if (!skip(myName)) {
        processInputs = channel.empty()
        params.tads.each {
            planName, analysisPlan ->
            analysisPlan["plan_name"] = planName
            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.mcool, analysisPlan)
                }
                | concat(processInputs)
                | set{processInputs}
        }
        // Note: nf-test raises error if calling process within
        // params.loops.each. Don't know why. Collect results and call
        // outside of each closure.
        processInputs | TADS
    }
    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}