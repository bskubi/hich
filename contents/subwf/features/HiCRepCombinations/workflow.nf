include {createCompositeStrategy; filterSamplesByStrategy; groupSamplesByStrategy} from '../../util/analysisPlans.nf'
include {skip; emptyOnLastStep} from '../../util/cli.nf'
include {columns} from '../../util/reshape.nf'
include {HICREP_COMBINATIONS} from './process.nf'

workflow HiCRepCombinations {
    take:
    samples

    main:
    myName = "HiCRepCombinations"
    if (!skip(myName)) {
        processInputs = channel.empty()

        params.hicrep.each {
            planName, analysisPlan ->
            analysisPlan["plan_name"] = planName
            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)
            filtered = filterSamplesByStrategy(samples, strategy)            
            grouped = groupSamplesByStrategy(filtered, strategy)
            grouped
                | filter{it.size() >= 2}
                | map{columns(it, ["dropAllNull":true])}
                | map{
                    samples ->
                    tuple(planName, samples.mcool, analysisPlan)
                }
                | concat(processInputs)
                | set{processInputs}
        }

        processInputs | HICREP_COMBINATIONS
    }

    emit:
    samples
}