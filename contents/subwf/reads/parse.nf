include {QCPairs} from './qcPairs.nf'
include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {PairtoolsParse2} from './parse/pairtoolsParse2.nf'


workflow Parse {
    take:
    samples

    main:
    if (!skip("parse")) {
        samples
            | filter{it.datatype in ["fastq", "sambam"]}
            | map{tuple(it.id, it.sambam, it.chromsizes, it.assembly, it.pairtoolsParse2Params, it.parseSQL, it.minMapq)}
            | PairtoolsParse2
            | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
            | set{result}
        keyUpdate(samples, result, "id") | set{samples}

        // TODO: simplify these workflow control steps since they
        // are repeated frequently.
        if ("parse" in params.general.get("qcAfter")) {
            QCPairs(samples, ["pairs"], "parse")
        }
    }

    samples = emptyOnLastStep("parse", samples)

    emit:
        samples
}
