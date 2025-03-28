include {groupRowsToColumnFormat; coalesce} from '../util/group.nf'
include {makeID} from '../util/samples.nf'
include {pack} from '../util/join.nf'


def processMerge (samplesToMerge, groupAttributes, proc, level, samplesAlreadyMerged) {
    def groupsToMerge = groupRowsToColumnFormat(samplesToMerge, groupAttributes, ["dropNull":true])
    
    def merged = groupsToMerge
        | map{coalesce(it)}
        | map{it += [id: makeID(it, columns = true)]}
        | map{tuple(it.id, it.latestPairs)}
        | proc
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}

    def attributes = groupsToMerge
        | map{coalesce(it, "_drop")}
        | map{it += [id: makeID(it, columns = true), aggregateLevel: level]}
   
    pack(attributes, merged)
    | concat(samplesAlreadyMerged)
}