
include {LabelAggregationPlans} from './LabelAggregationPlans/workflow.nf'
include {MergeTechrepsToBioreps} from './MergeTechrepsToBioreps/workflow.nf'
include {DedupPairs} from './DedupPairs/workflow.nf'
include {MergeBiorepsToConditions} from './MergeBiorepsToConditions/workflow.nf'
include {SplitPairs} from './SplitPairs/workflow.nf'
include {QCPairs} from '../reads/QCPairs/workflow.nf'
include {emptyOnLastStep; skip} from '../util/cli.nf'


workflow AggregatePairs {
    take:
    samples

    main:
    myName = "AggregatePairs"

    if (!skip(myName)) {
        samples
            | LabelAggregationPlans
            | MergeTechrepsToBioreps
            | DedupPairs
            | MergeBiorepsToConditions
            | SplitPairs
            | set{samples}

        if (myName in params.general.get("qcAfter")) {
            QCPairs(samples, ["latestPairs"], myName)
        }
    }

    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}