include {MergePairs} from '../MergePairs/workflow.nf'
include {LabelAggregationPlans} from '../LabelAggregationPlans/workflow.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'

workflow MergeBiorepsToConditions {
    take:
    samples

    main:
    myName = "MergeBiorepsToConditions"
    if (!skip(myName) && !skip("MergePairs")) {
        samples
            | branch {
                yes: it.aggregateLevel == "biorep" && !it.skipMerge && !it.skipBiorepMerge && it.mergeBiorepToCondition
                no: true
            }
            | set {mergeBioreps}

        MergePairs(
            mergeBioreps.yes, 
            ["cell", "condition", "aggregationPlanName"],
            "condition"
            )
            | LabelAggregationPlans
            | concat(samples)
            | set{samples}
    }

    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}