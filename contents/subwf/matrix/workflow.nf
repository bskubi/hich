include {LabelMatrixPlans} from './LabelMatrixPlans/workflow.nf'
include {HicMatrix} from './HicMatrix/workflow.nf'
include {McoolMatrix} from './McoolMatrix/workflow.nf'
include {IngestMatrix} from './IngestMatrix/workflow.nf'
include {emptyOnLastStep; skip} from '../util/cli.nf'


workflow CreateMatrix {
    take:
    samples

    main:
    myName = "CreateMatrix"

    if (!skip(myName)) {
        samples
            | LabelMatrixPlans
            | HicMatrix
            | McoolMatrix
            | IngestMatrix
            | set{samples}
    }

    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}