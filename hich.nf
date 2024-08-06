include {LoadSampleFile} from './subwf/loadSampleFile.nf'
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
    // we need to give a conda-based option for all workflow steps if possible
    // add read downsample step after select (can also be used for ingestion)

    // is there a way to clean up SLURM output from Nextflow?

    // compute the md5hash for the samples.csv and nextflow.config and save
    // as that hash in something like "runs" as a record of the analysis
    // actually it would be better to save a snapshot of params when nextflow
    // is launched, since its value is the result of aggregating several
    // places where param values can be set.

    // ingest .mcool/.hic files and interconvert between them using hictk

      

    LoadSampleFile
        | AssignParams
        | HeadReads
        | Align
        | Parse
        | IngestPairs
        | OptionalFragtag
        // | TechrepsToBioreps
        | Deduplicate
    //     | BiorepsToConditions
        | Select
    //     | MakeHic
        | MakeMcool
        | Hicrep
    //     | CallCompartments
    //     | CallLoops
    //     | CallInsulation
        
}

