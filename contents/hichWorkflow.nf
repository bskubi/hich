include {ParseParams} from './subwf/setup/parseParams.nf'
include {Setup} from './subwf/setup/setup.nf'

include {Align} from './subwf/reads/Align/workflow.nf'
include {ParseToPairs} from './subwf/reads/ParseToPairs/workflow.nf'
include {TagRestrictionFragments} from './subwf/reads/TagRestrictionFragments/workflow.nf'
include {IngestPairs} from './subwf/reads/IngestPairs/workflow.nf'
include {SelectPairs} from './subwf/reads/SelectPairs/workflow.nf'

include {AggregatePairs} from './subwf/aggregate/workflow.nf'

include {CreateMatrix} from './subwf/matrix/workflow.nf'
include {HiCRepCombinations} from './subwf/features/HiCRepCombinations/workflow.nf'
include {CompartmentScores} from './subwf/features/CompartmentScores/workflow.nf'
include {InsulationScores} from './subwf/features/InsulationScores/workflow.nf'
include {TADs} from './subwf/features/TADs/workflow.nf'
include {Loops} from './subwf/features/Loops/workflow.nf'
include {DifferentialLoops} from './subwf/features/DifferentialLoops/workflow.nf'

include {emptyOnLastStep; skip} from './subwf/util/cli.nf'

workflow HichWorkflow {
    
    ParseParams
        | Setup

        | Align
        | ParseToPairs
        | IngestPairs
        | TagRestrictionFragments
        | SelectPairs

        | AggregatePairs
        | CreateMatrix

        | HiCRepCombinations
        | CompartmentScores
        | InsulationScores
        | TADs
        | Loops
        | DifferentialLoops

        | set{samples}

    emit:
    samples
}

