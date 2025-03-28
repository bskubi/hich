include {QCReads} from './qcHicReads.nf'
include {withLog; stubLog; emptyOnLastStep; pack; isExistingFile; skip} from '../extraops.nf'

process HichFragtag {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null},
               mode: params.general.publish.mode
    container params.general.hichContainer
    label 'doJobArray'
    label 'pairs'

    input:
    tuple val(id), path(pairs), path(fragmentIndex)

    output:
    tuple val(id), path("${id}_fragtag.pairs.gz")

    shell:
    cmd = "hich fragtag '${fragmentIndex}' '${id}_fragtag.pairs.gz' '${pairs}'"

    logMap = [
        task: "HichFragtag",
        input: [
            id: id,
            pairs: pairs,
            fragmentIndex: fragmentIndex
        ],
        output: [
            pairs: "${id}_fragtag.pairs.gz"
        ]
    ]

    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}_fragtag.pairs.gz'"

    cmd = "hich pairs bin --idx1 rfrag1 --start1 rfrag_start1 --end1 rfrag_end1 --idx2 rfrag2 --start2 rfrag_end2 --end2 rfrag_end2 '${pairs}' '${fragmentIndex}' '${id}_fragtag.pairs.gz'"

    logMap = [
        task: "HichFragtag",
        input: [
            id: id,
            pairs: pairs,
            fragmentIndex: fragmentIndex
        ],
        output: [
            pairs: "${id}_fragtag.pairs.gz"
        ]
    ]

    stubLog(stub, cmd, logMap)
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