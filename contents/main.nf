include {ParseParams} from './subwf/setup/parseParams.nf'
include {Setup} from './subwf/setup/setup.nf'

include {Align} from './subwf/reads/align.nf'
include {FastqHead} from './subwf/reads/fastqHead.nf'
include {Parse} from './subwf/reads/parse.nf'
include {TagRestrictionFragments} from './subwf/reads/tagRestrictionFragments.nf'
include {IngestPairs} from './subwf/reads/ingestPairs.nf'
include {Select} from './subwf/reads/select.nf'

include {LabelAggregationPlans} from './subwf/aggregate/labelAggregationPlans.nf'
include {Aggregate} from './subwf/aggregate/aggregate.nf'

include {HicMatrix} from './subwf/matrix/hicMatrix.nf'
include {McoolMatrix} from './subwf/matrix/mcoolMatrix.nf'
include {IngestMatrix} from './subwf/matrix/ingestMatrix.nf'
include {Hicrep} from './subwf/features/hicrep.nf'
include {CompartmentScores} from './subwf/features/compartmentScores.nf'
include {Loops} from './subwf/features/loops.nf'
include {DifferentialLoops} from './subwf/features/differentialLoops.nf'
include {InsulationScores} from './subwf/features/insulationScores.nf'
include {emptyOnLastStep; skip} from './subwf/util/cli.nf'

workflow {
    
    ParseParams              // Setup workflow inputs
        | Setup
        | FastqHead

        | Align                 // Align .fastq -> .bam
        | Parse                 // .bam -> .pairs and read-level filters
        | IngestPairs
        | TagRestrictionFragments   // Only affects samples with "fragmentIndex"
        | Select

        | LabelAggregationPlans
        | Aggregate

        | HicMatrix             // Create contact matrices
        | McoolMatrix
        | IngestMatrix

        | Hicrep                // Call features
        | CompartmentScores
        | Loops
        | DifferentialLoops
        | InsulationScores

        | set{samples}

    samples = emptyOnLastStep("end", samples)
}

