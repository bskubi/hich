include {QCPairs} from '../qcPairs.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {PARSE_TO_PAIRS} from './process.nf'


workflow ParseToPairs {
    take:
    samples

    main:
    if (!skip("ParseToPairs")) {
        samples
            | filter{it.datatype in ["fastq", "sambam"]}
            | map{tuple(it.id, it.sambam, it.chromsizes, it.assembly, it.pairtoolsParse2Params, it.parseSQL, it.minMapq)}
            | PARSE_TO_PAIRS
            | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
            | set{result}
        keyUpdate(samples, result, "id") | set{samples}

        // TODO: simplify these workflow control steps since they are repeated frequently.
        if ("ParseToPairs" in params.general.get("qcAfter")) {
            QCPairs(samples, ["pairs"], "ParseToPairs")
        }
    }

    samples = emptyOnLastStep("ParseToPairs", samples)

    emit:
        samples
}
