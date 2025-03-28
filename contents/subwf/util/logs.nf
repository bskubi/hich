import groovy.json.JsonOutput

def prepForJson(obj) {
    if (obj instanceof Map) {
        return obj.collectEntries { k, v ->
            [(k.toString()): prepForJson(v)]
        }
    } else if (obj.getClass() in [List, ArrayList]) {
        return obj.withIndex().collectEntries { item, idx ->
            [(idx): prepForJson(item)]
        }
    } else {
        return obj.toString()
    }
}


def withLog(cmdArg, mapArg, stubArg = null) {
    def map = mapArg + [command: cmdArg.toString()]
    if (stubArg) {
        map += [stub: stubArg.toString()]
    }
    def useArg = stubArg ? stubArg : cmdArg

    // Convert to a map to avoid cyclic entries
    def preppedMap = prepForJson(map)

    def json = JsonOutput.toJson(preppedMap)

    return "${useArg} && cat <<'EOF' > hich_output.json\n${json}\nEOF"
}

def stubLog(stubArg, cmdArg, mapArg) {
    def map = mapArg + [command: cmdArg.toString(), stub: stubArg.toString()]
    // Convert to a map to avoid cyclic entries
    def preppedMap = prepForJson(map)
    def json = JsonOutput.toJson(preppedMap)
    return "${stubArg} && cat <<'EOF' > hich_output.json\n${json}\nEOF"
}