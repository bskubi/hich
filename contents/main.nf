include {ParseArgs} from './subwf/parseArgs.nf'
include {Setup} from './subwf/setup.nf'
include {Align} from './subwf/align.nf'
include {FastqHead} from './subwf/fastqHead.nf'
include {Parse} from './subwf/parse.nf'
include {TagRestrictionFragments} from './subwf/tagRestrictionFragments.nf'
include {IngestPairs} from './subwf/ingestPairs.nf'
include {Select} from './subwf/select.nf'
include {HicMatrix} from './subwf/hicMatrix.nf'
include {McoolMatrix} from './subwf/mcoolMatrix.nf'
include {IngestMatrix} from './subwf/ingestMatrix.nf'
include {Hicrep} from './subwf/hicrep.nf'
include {CompartmentScores} from './subwf/compartmentScores.nf'
include {Loops} from './subwf/loops.nf'
include {DifferentialLoops} from './subwf/differentialLoops.nf'
include {InsulationScores} from './subwf/insulationScores.nf'
include {emptyOnLastStep; skip} from './subwf/extraops.nf'
include {AggregateTechreps; AggregateBioreps; AggregateConditions} from './subwf/aggregate.nf'

workflow {

    ParseArgs              // Setup workflow inputs
        | Setup
        | FastqHead

        | Align                 // Align .fastq -> .bam
        | Parse                 // .bam -> .pairs and read-level filters
        | IngestPairs
        | TagRestrictionFragments   // Only affects samples with "fragmentIndex"
        | Select

        | AggregateTechreps     // Downsample, deduplicate, merge
        | AggregateBioreps
        | AggregateConditions

        | HicMatrix             // Create contact matrices
        | McoolMatrix
        | IngestMatrix

        | Hicrep                // Call features
        | CompartmentScores
        | Loops
        | DifferentialLoops
        | InsulationScores

        | set{samples}

    samples = emptyOnLastStep("End", samples)
}

