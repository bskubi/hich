include {QCPairs} from '../QCPairs/workflow.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {SELECT_PAIRS} from './process.nf'

workflow SelectPairs {
    take:
    samples
    
    main:

    if (!skip("Select")) {
        samples
            | filter{it.latestPairs && it.selectPairs_opts && it.selectPairs_opts?.filters}
            | map{tuple(it.id, it.latestPairs, it.selectPairs_opts)}
            | SELECT_PAIRS
            | map{[id:it[0], selectPairs:it[1], latest:it[1], latestPairs:it[1]]}
            | set{result}
        keyUpdate(samples, result, "id") | set{samples}

        if ("SelectPairs" in params.general.get("qcAfter")) {
            QCPairs(samples, ["selectPairs"], "SelectPairs")
        }
    }

    samples = emptyOnLastStep("SelectPairs", samples)

    emit:
    samples
}