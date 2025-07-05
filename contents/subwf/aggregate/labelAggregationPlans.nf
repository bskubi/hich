include {columnsToRows} from '../util/reshape.nf'
include {emptyOnLastStep} from '../util/cli.nf'
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

        // Add aggregation plans to samples not yet labeled with one
        channel.of(params.aggregate)
            | map {[
                    planName: it.keySet() as List, 
                    planParams: it.values() as List
                ]}
            | columnsToRows
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