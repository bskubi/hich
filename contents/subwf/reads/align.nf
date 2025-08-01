include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {withLog; stubLog} from '../util/logs.nf'
include {AlignBwa} from './processes/alignBwa'
include {getFastq} from './helpers/alignHelpers.nf'

workflow Align {
    take:
    samples

    main:

    if (!skip("align")) {

        samples
            | filter {it.datatype == "fastq"}
            | map{tuple(it.id, it.aligner, it.alignerIndexDir, it.alignerIndexPrefix, getFastq(it), it.bwaFlags, it.minMapq)}
            | AlignBwa
            | map{[id:it[0], sambam:it[1], latest:it[1], latestSambam:it[1]]}
            | set{results}
    

        keyUpdate(samples, results, "id") | set{samples}
    }


    samples = emptyOnLastStep("align", samples)

    emit:
    samples
}

