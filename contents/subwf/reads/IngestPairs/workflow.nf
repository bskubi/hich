include {QCPairs} from '../QCPairs/workflow.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {INGEST_PAIRS} from './process.nf'

workflow IngestPairs {
    take:
    samples

    main:

    if (!skip("IngestPairs")) {
        samples
            | filter{it.datatype == "pairs"}
            | map{tuple(it.id, it.pairs, it.chromsizes, it.reshapeParams)}
            | INGEST_PAIRS
            | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
            | set{result}

        keyUpdate(samples, result, "id") | set{samples}

        if ("IngestPairs" in params.general.get("qcAfter")) {
            samples = QCPairs(samples, ["pairs"], "IngestPairs")
        }
    }

    samples = emptyOnLastStep("IngestPairs", samples)

    emit:
    samples
}