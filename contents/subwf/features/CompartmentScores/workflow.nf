include {createCompositeStrategy; filterSamplesByStrategy} from '../../util/analysisPlans.nf'
include {skip; emptyOnLastStep} from '../../util/cli.nf'
include {COMPARTMENT_SCORES} from './process.nf'

workflow CompartmentScores {
    take:
    samples

    main:
    
    myName = "CompartmentScores"
    if (!skip(myName)) {
        processInputs = channel.empty()

        params.compartments.each {
            planName, analysisPlan ->
            analysisPlan["plan_name"] = planName
            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)
            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.genomeReference, sample.mcool, analysisPlan)
                }
                | concat(processInputs)
                | set{processInputs}
        }
        // Note: nf-test raises error if calling process within
        // params.loops.each. Don't know why. Collect results and call
        // outside of each closure.
        processInputs | COMPARTMENT_SCORES
    }
    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}
