include {createCompositeStrategy; filterSamplesByStrategy} from '../../util/analysisPlans.nf'
include {skip; emptyOnLastStep} from '../../util/cli.nf'
include {INSULATION_SCORES} from './process.nf'
workflow InsulationScores {
    take:
    samples

    main:

    myName = "InsulationScores"
    if (!skip(myName)) {
        processInputs = channel.empty()
        params.insulation.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.mcool, analysisPlan.resolution, analysisPlan.cooltoolsInsulationParams, analysisPlan.window)
                }
                | concat(processInputs)
                | set{processInputs}
        }
        processInputs | INSULATION_SCORES
    }
    samples = emptyOnLastStep(myName, samples)


    emit:
    samples
}