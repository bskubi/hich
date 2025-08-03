include {GroupToColumnar} from '../../util/groupToColumnar.nf'
include {MERGE_PAIRS} from './process.nf'
include {makeID} from '../../util/samples.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {coalesce} from '../../util/reshape.nf'

workflow MergePairs {
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
    )
        | map{
            if (level == "biorep") {
                result = it.findAll {key, value -> !(key in ["techrep", "id"])}
            }
            else if (level == "condition") {
                result = it.findAll {key, value -> !(key in ["techrep", "biorep", "id"])}
            }
            result
        }
        | set{sampleGroups}
   
    // Create new id and pairs attributes for merged sample
    sampleGroups
        | map{coalesce(it, false)}
        | map{it += [id: makeID(it, true)]}
        | map{tuple(it.id, it.latestPairs instanceof List ? it.latestPairs : [it.latestPairs])}
        | MERGE_PAIRS
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{merged}

    // Obtain attribute values for the merged sample 
    sampleGroups
        | map{coalesce(it, true)}   // Common attributes from input samples
        | map{it += [id: makeID(it, columns = false), aggregateLevel: level]} // New id, aggregateLevel
        | set{attributes}

    keyUpdate(attributes, merged, "id") | set{result}

    emit:
    result
}
