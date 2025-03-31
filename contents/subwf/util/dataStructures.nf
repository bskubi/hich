def asHashSet(val) {
    // Convert a non-HashSet val into a HashSet
    if (val instanceof ArrayList) return new HashSet(val)
    if (val instanceof HashSet) return val
    return new HashSet([val])
}

def asList(val) {
    return (val instanceof ArrayList ? val : [val])
}

def cmdList(val, delimiter = " ", quotes = true) {
    delimiter = quotes ? "\"${delimiter}\"" : delimiter
    def joined = asList(val).join(delimiter)
    return quotes ? "\"${joined}\"" : joined
}