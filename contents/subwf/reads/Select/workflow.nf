include {QCPairs} from '../qcPairs.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {SELECT} from './process.nf'

workflow Select {
    take:
    samples
    
    main:

    if (!skip("Select")) {
        samples
            | filter{(it.pairtoolsSelectParams || it.pairtoolsSelectFilters) && (it.pairs || it.latestPairs)}
            | map{tuple(it.id, it.latestPairs, it.pairtoolsSelectParams, it.pairtoolsSelectFilters)}
            | SELECT
            | map{[id:it[0], selectPairs:it[1], latest:it[1], latestPairs:it[1]]}
            | set{result}
        keyUpdate(samples, result, "id") | set{samples}

        if ("Select" in params.general.get("qcAfter")) {
            QCPairs(samples, ["selectPairs"], "Select")
        }
    }

    samples = emptyOnLastStep("select", samples)

    emit:
    samples
}