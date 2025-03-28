import groovy.json.JsonOutput

def prepForJson(obj) {
    def updated = obj
    if (obj instanceof Map) {
        updated = [:]
        obj.each {
            k, v ->
            kStr = k.toString()
            updated[(kStr)] = prepForJson(v)
        }
    } else if (obj.getClass() in [List, ArrayList]) {
        updated = []
        obj.eachWithIndex {it, idx -> updated[idx] = prepForJson(it) }
    } else {
        updated = obj.toString()
    }
    return updated
}


def withLog(cmdArg, mapArg) {
    def map = mapArg + [command: cmdArg.toString()]
    // Convert to a map to avoid cyclic entries
    def preppedMap = prepForJson(map)
    def json = JsonOutput.toJson(preppedMap)
    return "${cmdArg} && cat <<'EOF' > hich_output.json\n${json}\nEOF"
}

def stubLog(stubArg, cmdArg, mapArg) {
    def map = mapArg + [command: cmdArg.toString(), stub: stubArg.toString()]
    // Convert to a map to avoid cyclic entries
    def preppedMap = prepForJson(map)
    def json = JsonOutput.toJson(preppedMap)
    return "${stubArg} && cat <<'EOF' > hich_output.json\n${json}\nEOF"
}