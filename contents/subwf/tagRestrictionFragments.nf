include {QCReads} from './qcHicReads.nf'
include {emptyOnLastStep; pack; isExistingFile; skip} from './extraops.nf'

process HichFragtag {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"
    label 'doJobArray'
    label 'pairs'

    input:
    tuple val(id), path(pairs), path(fragmentIndex)

    output:
    tuple val(id), path("${id}_fragtag.pairs.gz")

    shell:
    "hich fragtag ${fragmentIndex} ${id}_fragtag.pairs.gz ${pairs}"

    stub:
    "touch ${id}_fragtag.pairs.gz"
}

workflow TagRestrictionFragments {
    take:
        samples

    main:

    samples
        | filter{!skip("tagRestrictionFragments") && isExistingFile(it.fragmentIndex) && isExistingFile(it.latestPairs)}
        | map{tuple(it.id, it.latestPairs, it.fragmentIndex)}
        | HichFragtag
        | map{[id:it[0], fragPairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{result}
    pack(samples, result) | set{samples}

    if ("tagRestrictionFragments" in params.general.get("qcAfter")) {
        samples = QCReads(samples, "tagRestrictionFragments")
    }

    samples = emptyOnLastStep("tagRestrictionFragments", samples)

    emit:
        samples
}