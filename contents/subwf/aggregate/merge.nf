include {GroupToColumnar} from '../util/groupToColumnar.nf'
include {PairtoolsMerge} from './processes/pairtoolsMerge.nf'
include {makeID} from '../util/samples.nf'
include {keyJoin} from '../util/keyJoin.nf'
include {coalesce} from '../util/reshape.nf'

workflow Merge {
    take:
    samples
    groupBy
    level

    main:

    // Group samples by groupBy keys (i.e. condition, biorep, techrep) and convert
    // to columnar format after values by id to get deterministic output so -resume will work. 
    GroupToColumnar(
        samples, 
        groupBy, 
        ["id"], 
        ["dropAllNull":true]
    ) | set{sampleGroups}
    
    // Create new id and pairs attributes for merged sample
    sampleGroups
        | map{coalesce(it, false)}
        | map{it += [id: makeID(it, true)]}
        | map{tuple(it.id, it.latestPairs)}
        | PairtoolsMerge
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{merged}

    // Obtain attribute values for the merged sample 
    sampleGroups
        | map{coalesce(it, true)}   // Common attributes from input samples
        | map{it += [id: makeID(it, columns = false), aggregateLevel: level]} // New id, aggregateLevel
        | set{attributes}

    keyJoin(attributes, merged, "id") | set{result}

    emit:
    result
}