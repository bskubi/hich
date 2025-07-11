include {QCReads} from './qcHicReads.nf'
include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {isExistingFile} from '../util/files.nf'
include {withLog; stubLog} from '../util/logs.nf'

process HichFragtag {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    conda "$projectDir/env/dev_env.yml"

    input:
    tuple val(id), path(pairs), path(fragmentIndex)

    output:
    tuple val(id), path(tagged)

    shell:
    tagged = "${id}_fragtag.pairs.gz"
    cmd = "hich pairs map-ends --idx1 rfrag1 --start1 rfrag_start1 --end1 rfrag_end1 --idx2 rfrag2 --start2 rfrag_start2 --end2 rfrag_end2 '${fragmentIndex}' '${pairs}' '${tagged}'"

    logMap = [
        task: "HichFragtag",
        input: [
            id: id,
            pairs: pairs,
            fragmentIndex: fragmentIndex
        ],
        output: [
            pairs: tagged
        ]
    ]

    withLog(cmd, logMap)

    stub:
    tagged = "${id}_fragtag.pairs.gz"
    stub = "touch '${tagged}'"

    cmd = "hich pairs map-ends --idx1 rfrag1 --start1 rfrag_start1 --end1 rfrag_end1 --idx2 rfrag2 --start2 rfrag_start2 --end2 rfrag_end2 '${fragmentIndex}' '${pairs}' '${tagged}'"

    logMap = [
        task: "HichFragtag",
        input: [
            id: id,
            pairs: pairs,
            fragmentIndex: fragmentIndex
        ],
        output: [
            pairs: tagged
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
    keyUpdate(samples, result, "id") | set{samples}

    if ("tagRestrictionFragments" in params.general.get("qcAfter")) {
        samples = QCReads(samples, "tagRestrictionFragments")
    }

    samples = emptyOnLastStep("tagRestrictionFragments", samples)

    emit:
        samples
}