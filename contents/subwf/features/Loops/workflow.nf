include {createCompositeStrategy; filterSamplesByStrategy;} from '../../util/analysisPlans.nf'
include {skip; emptyOnLastStep} from '../../util/cli.nf'
include {LOOPS} from './process.nf'

workflow Loops {
    take:
    samples

    main:
    myName = "Loops"
    if (!skip(myName)) {
        processInputs = channel.empty()

        params.loops.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.latestMatrix, analysisPlan.mustacheParams)
                }
                | set{processInputs}
        }
        // Note: nf-test raises error if calling process within
        // params.loops.each. Don't know why. Collect results and call
        // outside of each closure.
        processInputs | LOOPS
    }

    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}