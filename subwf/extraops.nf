include {keydiff; sqljoin} from './sqljoin.nf'

def toHashMap(keys, vals) {
    if (!(vals instanceof ArrayList) || vals.size() == 1) {
        return [keys:vals]
    }
    return [keys, vals].transpose()
                       .collectEntries { [it[0], it[1]] }
}

workflow transact {
    take:
        ch
        proc
        input
        output
        tags

    main:
        ch
            | map {
                hashmap ->
                proc_inputs = input.collect{key -> hashmap.get(key)}
                proc_inputs.size() == 1 ? proc_inputs[0] : tuple(*proc_inputs)
            }
            | proc
            | map{
                proc_outputs ->
                result_map = toHashMap(output, proc_outputs)
                tags.each {
                    tag, orig ->
                    result_map[tag] = result_map[orig]
                }
                
                result_map
            }
            | set{result}
    emit:
    result
}

def getIdx(item, index) {
    /* If item is a list, get value at index, where too-large indexes are set
       to index of final item. If item is not a list, return it.
    */
    if (item instanceof List) {
        def validIndex = Math.min(index, item.size() - 1)
        return item[validIndex]
    }
    return item
}

workflow pack {
    take:
    channels
    join_by

    main:
    
    result = null
    skip = channels[1..-1]
    skip.eachWithIndex {
        cur_ch, i ->

        prev_ch = channels[i]
        next_ch = channels[i+1]

        by = getIdx(join_by, i)

        skip[i] = sqljoin(cur_ch,
                          prev_ch,
                          [by:by, suffix:""])
        result = skip[i]
    }

    emit:
    result
}

workflow transpack {
    take:
    ch
    proc
    input
    output
    tags
    chPack
    joinBy

    main:
    transact(ch, proc, input, output, tags) | set {obtained}
    chPack = chPack instanceof List ? chPack : [chPack]
    toPack = [obtained] + chPack
    pack(toPack, joinBy) | set{result}

    emit:
    result
}