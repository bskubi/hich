include {keydiff; sqljoin} from './sqljoin.nf'

def toHashMap(keys, vals) {
    if (!(vals instanceof ArrayList) || vals.size() == 1) {
        return [keys:vals]
    }
    return [keys, vals].transpose()
                       .collectEntries { [it[0], it[1]] }
}

def transact (proc, ch, input, output, tags = [:]) {
    /* Extract process inputs from hashmap channel items, call the process,
    and rebuild a hashmap with process outputs. "Tag" some outputs by assigning
    them to a second key.
    */
    return ch
            | map {
                hashmap ->
                def proc_inputs = input.collect{key -> hashmap.get(key)}
                proc_inputs.size() == 1 ? proc_inputs[0] : tuple(*proc_inputs)
            }
            | proc
            | map{
                proc_outputs ->
                def result_map = toHashMap(output, proc_outputs)
                tags.each {
                    tag, orig ->
                    result_map[tag] = result_map[orig]
                }
                
                result_map
            }
}

def pack(channels, joinBy = "id") {
    // Extract first and remaining channels
    def first = channels[0]
    def rest = channels[1..-1]

    // Iteratively join the channels from first to last.
    // If joinBy is a list, join on the joinBy element corresponding to each
    // item. Otherwise, join on joinBy every time.
    sizeMatch = !(joinBy instanceof List)
                || joinBy.size() == channels.size() - 1
    assert sizeMatch, "In extraops.nf, there must be a single non-list joinBy or one joinBy for every join"
    return rest.inject(
        first,
        {
            addMe, addTo ->
            // Get the key to join on
            idxOfAddMe = channels.indexOf(addMe)
            def by = joinBy instanceof List ? joinBy[idxOfAddMe] : joinBy

            // Join the previous results (addMe) to the new (presumably larger)
            // hashmaps in addTo. On overlapping columns for a given
            // channel item, addMe's values replace addTo's values.
            sqljoin(addTo, addMe, [by: by, suffix: ""])
        }
    )
}

def transpack (proc, channels, input, output, tags = [:], by = "id") {
    // Convenience function to call transact followed by pack.
    def channels_list = channels instanceof List ? channels : [channels]
    def obtained = transact(proc, channels_list[0], input, output, tags)
    return pack([obtained] + channels_list, by)
}

import java.nio.file.Path
import java.nio.file.Paths


def hashmapdiff(ch1, ch2, by, how = "left", suffix = "__joindiff__") {
    def joined = sqljoin(ch1, ch2, [by: by, how:how, suffix:"__joindiff__"])
    def keyIsInBy = {k -> by instanceof List ? k in by : k == by}

    return joined.collect {
        hashmap ->

        def diff = [:]
        hashmap.each {
            k1, v1 ->
            if (!k1.endsWith("__joindiff__") && !keyIsInBy(k1)) {
                def k2 = k1 + "__joindiff__"
                if (!k2 in hashmap.keySet()) {
                    diff[k1] = [v1]
                } else if (hashmap[k2] != v1) {
                    v2 = hashmap[k2]
                    if (v1 instanceof Path && v2 instanceof Path) {
                        if (v1.getFileName().toString() != v2.getFileName().toString()) {
                            diff[k1] = [v1, hashmap[k2]]
                        }
                    } else {
                        diff[k1] = [v1, hashmap[k2]]
                    }
                    
                }
            }
        }
        diff
    }
}