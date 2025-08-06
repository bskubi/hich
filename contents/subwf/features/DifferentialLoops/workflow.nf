include {createCompositeStrategy; pairSamplesByStrategy} from '../../util/analysisPlans.nf'
include {skip; emptyOnLastStep} from '../../util/cli.nf'
include {DIFFERENTIAL_LOOPS} from './process.nf'

workflow DifferentialLoops {
    take:
    samples

    main:
    myName = "DifferentialLoops"
    if (!skip(myName)) {
        processInputs = channel.empty()

        params.diffloops.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)
            pairSamplesByStrategy(samples, strategy)
                | filter{s1, s2 -> s1 != s2 && s1.matrixPlanName == s2.matrixPlanName}
                | map{
                    s1, s2 ->
                    prefix = "${s1.id}_${s2.id}"
                    tuple(s1.id, s2.id, prefix, s1.latestMatrix, s2.latestMatrix, analysisPlan)
                }
                | concat(processInputs)
                | set{processInputs}
        }
        processInputs | DIFFERENTIAL_LOOPS
    }
    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}