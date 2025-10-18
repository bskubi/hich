include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {withLog; stubLog} from '../../util/logs.nf'
include {ALIGN} from './process.nf'
include {getFastq} from './functions.nf'
include {extract_inputs; load_sample_attrib_schema} from '../../util/schema.nf'

workflow Align {
    take:
    samples

    main:
    schema = load_sample_attrib_schema()

    if (!skip("Align")) {

        samples
            | filter {it.datatype == "fastq"}
            | map{
                sample ->
                inputs = extract_inputs("ALIGN", schema, sample)
                inputs += ["fastq": inputs.subMap("fastq1", "fastq2", "fastq").values()]
                inputs.subMap("id", "aligner", "aligner_index_dir", "aligner_index_prefix", "fastq", "align_opts", "min_mapq").values().toList()
            }
            | ALIGN
            | map{[id:it[0], sambam:it[1], latest:it[1], latestSambam:it[1]]}
            | set{results}
        keyUpdate(samples, results, "id") | set{samples}
    }


    samples = emptyOnLastStep("Align", samples)

    emit:
    samples
}

