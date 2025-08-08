include {DEDUP_PAIRS} from './process.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {QCPairs} from '../../reads/QCPairs/workflow.nf'

workflow DedupPairs {
    take:
    samples

    main:
    myName = "DedupPairs"

    // Deduplicate techreps, bioreps, and input conditions
    if (!skip(myName)) {
        samples
        | branch {
            yes: (
                (it.aggregateLevel in ["techrep", "biorep"])
                && !(it.aggregateLevel == "techrep" && it.skipDedupTechreps)
                && !(it.aggregateLevel == "biorep" && it.skipDedupBioreps)
                && !(it.aggregateLevel in ["techrep", "biorep"] && it.skipDedup)
            )
            no: true
        }
        | set{dedup}

        dedup.yes
        | map{tuple(it.id, it.latestPairs, it.dedupSingleCell, it.dedupMaxMismatch, it.dedupMethod, it.pairtoolsDedupParams)}
        | DEDUP_PAIRS
        | map{[id:it[0], dedupPairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set {deduplicated}

        keyUpdate(dedup.yes, deduplicated, "id")
        | concat(dedup.no)
        | set {samples}

        if (myName in params.general.get("qcAfter")) {
            QCPairs(samples, ["dedupPairs"], myName)
        }
    }

    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}