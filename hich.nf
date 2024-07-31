include {AssignParams} from './subwf/assignParams.nf'
include {Align} from './subwf/align.nf'
include {HeadReads} from './subwf/downsampleFastq.nf'
include {Parse} from './subwf/parseToContacts.nf'
include {OptionalFragtag} from './subwf/fragtag.nf'
include {TechrepsToBioreps; BiorepsToConditions} from './subwf/mergePairs.nf'
include {IngestPairs} from './subwf/ingestPairs.nf'
include {Deduplicate} from './subwf/dedup.nf'
include {Select} from './subwf/select.nf'
include {MakeMcool} from './subwf/makeMcool.nf'
include {MakeHic} from './subwf/makeHic.nf'
include {CallLoops} from './subwf/callLoops.nf'
include {CallCompartments} from './subwf/callCompartments.nf'
include {CallInsulation} from './subwf/callInsulation.nf'
include {Hicrep} from './subwf/hicrep.nf'

workflow {
    // publishDir mode should be easier to switch from "move" to "copy"
    // we need to give a conda-based option for all workflow steps if possible
    // add read downsample step after select (can also be used for ingestion)
    // is there a way to clean up SLURM output from Nextflow?
    sampleCSV = params.general.sampleCSV
    channel.fromPath(sampleCSV.filename, checkIfExists: true)
        | splitCsv(header: true, sep: sampleCSV.sep)
        | map{it.id = it.sample_id; it}
        | AssignParams
        | HeadReads
        | Align
        | Parse
        | IngestPairs
        | OptionalFragtag
        | TechrepsToBioreps
        | Deduplicate
        | BiorepsToConditions
        | Select
        | MakeHic
        | MakeMcool
        | Hicrep
        | CallCompartments
        | CallLoops
        | CallInsulation
        
}

