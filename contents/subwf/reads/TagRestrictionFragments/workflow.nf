include {QCPairs} from '../QCPairs/workflow.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {isExistingFile} from '../../util/files.nf'
include {TAG_RESTRICTION_FRAGMENTS} from './process.nf'

workflow TagRestrictionFragments {
    take:
        samples

    main:

    if (!skip("TagRestrictionFragments")) {
        samples
            | filter{isExistingFile(it.fragmentIndex) && isExistingFile(it.latestPairs)}
            | map{tuple(it.id, it.latestPairs, it.fragmentIndex, it.tag_restriction_fragments_opts)}
            | TAG_RESTRICTION_FRAGMENTS
            | map{[id:it[0], fragPairs:it[1], latest:it[1], latestPairs:it[1]]}
            | set{result}
        keyUpdate(samples, result, "id") | set{samples}

        if ("TagRestrictionFragments" in params.general.get("qcAfter")) {
            samples = QCPairs(samples, ["fragPairs"], "TagRestrictionFragments")
        }
    }


    samples = emptyOnLastStep("TagRestrictionFragments", samples)

    emit:
        samples
}