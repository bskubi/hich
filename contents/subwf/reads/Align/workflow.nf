include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {withLog; stubLog} from '../../util/logs.nf'
include {ALIGN} from './process.nf'
include {getFastq} from './functions.nf'

workflow Align {
    take:
    samples

    main:

    if (!skip("Align")) {

        samples
            | filter {it.datatype == "fastq"}
            | map{
                if (!it.alignerIndexDir) {
                    error("alignerIndexDir is null for sample ${it}")
                }
                if (!getFastq(it)) {
                    error("datatype is 'fastq' but no fastq files found for sample ${it}")
                }
            }
            | map{tuple(it.id, it.aligner, it.alignerIndexDir, it.alignerIndexPrefix, getFastq(it), it.align_opts, it.minMapq)}
            | ALIGN
            | map{[id:it[0], sambam:it[1], latest:it[1], latestSambam:it[1]]}
            | set{results}
    

        keyUpdate(samples, results, "id") | set{samples}
    }


    samples = emptyOnLastStep("Align", samples)

    emit:
    samples
}

