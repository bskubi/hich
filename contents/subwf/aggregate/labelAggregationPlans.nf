include {rowHashmapToRowChannel} from '../util/rowsCols.nf'
include {emptyOnLastStep} from '../util/workflowControl.nf'
include {Setup} from '../setup/setup.nf'

workflow LabelAggregationPlans {
    take:
    samples

    main:

    // Add aggregationPlanName and plan parameters to previously unlabeled samples
    if (params.containsKey("aggregate")) {

        // Split samples into those with and without an aggregation plan name
        samples
            | branch {
                yes: it.aggregationPlanName != null
                no: true
            }
            | set{alreadyLabeled}

        // For samples without an aggregation plan name, add one with the appropriate parameters
        // using the "aggregationPlans" section
        rowHashmapToRowChannel(params.aggregate, "planName", "planParams")
            | combine(alreadyLabeled.no)
            | map {
                plan = it[0]
                sample = it[1]
                sample += [aggregationPlanName: plan.planName] + plan.planParams
                sample
            }
            | set{newlyLabeled}
        
        alreadyLabeled.yes
            | concat(newlyLabeled)
            | set {samples}

        if (params.containsKey("keepUnaggregated")) {
            samples
                | concat(alreadyLabeled.no)
                | set{samples}
        }
    }

    emit:
    samples
}