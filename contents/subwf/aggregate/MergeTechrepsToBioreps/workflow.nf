include {MergePairs} from '../MergePairs/workflow.nf'
include {LabelAggregationPlans} from '../LabelAggregationPlans/workflow.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'

workflow MergeTechrepsToBioreps {
    take:
    samples

    main:
    myName = "MergeTechrepsToBioreps"

    if (!skip(myName) && !skip("MergePairs")) {
        samples
            | branch {
                yes: it.aggregateLevel == "techrep" && !it.skipMerge && !it.skipTechrepMerge && it.mergeTechrepToBiorep
                no: true
            }
            | set {mergeTechreps}


        // Merge techreps to bioreps
        MergePairs(
            mergeTechreps.yes, 
            ["cell", "biorep", "condition", "aggregationPlanName"],
            "biorep"
            )
        | LabelAggregationPlans
        | concat(samples)
        | set {samples}

    }

    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}