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

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.genomeReference, sample.mcool, analysisPlan.resolution, analysisPlan.eigsCisParams, analysisPlan.nEigs)
                }
                | concat(processInputs)
                | set{processInputs}
        }
        processInputs | COMPARTMENT_SCORES
    }
    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}
