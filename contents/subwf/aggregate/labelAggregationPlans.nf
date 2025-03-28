include {rowHashmapToRowChannel} from '../util/rowsCols.nf'
include {emptyOnLastStep} from '../util/workflowControl.nf'
include {Setup} from '../setup/setup.nf'

workflow LabelAggregationPlans {
    take:
    samples

    main:

    // Add aggregationPlanName and plan parameters to previously unlabeled samples
    if (params.containsKey("aggregationPlans")) {

        // Split samples into those with and without an aggregation plan name
        samples
            | branch {
                noAggPlan: it.aggregationPlanName == null
                hasAggPlan: true
            }
            | set{samples}

        // For samples without an aggregation plan name, add one with the appropriate parameters
        // using the "aggregationPlans" section
        rowHashmapToRowChannel(params.aggregationPlans, "planName", "planParams")
            | combine(samples.noAggPlan)
            | map {
                plan = it[0]
                sample = it[1]
                sample += [aggregationPlanName: plan.planName] + plan.planParams
                sample
            }
            | set{samples}

        if (params.containsKey("keepUnaggregated")) {
            samples | concat(samples.hasAggPlan) | set{samples}
        }
    }

    samples = emptyOnLastStep("labelAggregationPlans", samples)

    emit:
    samples
}