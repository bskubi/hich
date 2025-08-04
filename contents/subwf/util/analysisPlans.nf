include {asHashSet} from "./dataStructures.nf"

def createCompositeStrategy(strategyNames, sampleSelectionStrategies, combineHow = [:]) {
    /* 
        Example:

        strategyNames = ["strategy1", "strategy2"]
        sampleSelectionStrategies = [
            strategy1: [attribute1: [1, 2], attribute2: [3]]
            strategy2: [attribute2: [4], attribute3: [5, 6]]
        ]

        returns:
        [attribute1: [1, 2], attribute2: [4], attribute3: [5, 6]]
    */
    def strategiesToCombine = sampleSelectionStrategies.subMap(strategyNames).values() ?: []
    def compositeStrategy = [:]
    strategiesToCombine.each {
        strategy ->

        strategy.each {
            attribute, values ->

            def oldVals = compositeStrategy.get(attribute, new HashSet())
            def newVals = asHashSet(values)
            def updatedVals = combineHow.get(attribute, {a, b -> b})(oldVals, newVals)
            compositeStrategy += [(attribute): updatedVals]
        }
    }

    return compositeStrategy
}

def filterSamplesByStrategy(samples, strategy) {
    if (!strategy) return samples

    def filtered = samples | filter {
        sample ->
        strategy.every {
            attribute, keepValues ->
            def ignoreAttribute = (attribute in ["same", "different"])
            def keep = (sample.attribute in keepValues)
            print([attribute, ignoreAttribute, keep])
            ignoreAttribute || keep
        }
    }

    filtered
        | collect
        | {
            assert it.size() > 0, "Error: In filterSamplesByStrategy with sample selection strategy ${strategy}, no samples matched this filter."
        }

    return filtered
}

def groupSamplesByStrategy(samples, strategy) {
    def groupBy = strategy.groupBy ?: []
    return samples
        | map{
            tuple( it.subMap(groupBy), it )
        }
        | groupTuple
        | map{it[1]}
}

def pairSamplesByStrategy(samples, strategy) {

    def same = strategy.same ?: []
    def different = strategy.different ?: []
    strategy = strategy.findAll{key, value -> !(key in ["same", "different", "groupBy"])}
    def sameAndDifferent = same.intersect(different)
    if (!sameAndDifferent.isEmpty()) {
        error("Conflict in pairSamplesbyStrategy: 'same' and 'different' include same values ${sameAndDifferent}.")
    }

    def filtered = filterSamplesByStrategy(samples, strategy)

    def combined = filtered
        | combine(filtered)
        | filter{it[0].id <= it[1].id}
        | unique
        | filter {
            s1, s2 ->

            def sameOK = same.every {key -> s1.get(key) == s2.get(key)}
            def differentOK = different.every {key -> s1.get(key) != s2.get(key)}
            sameOK && differentOK
        }

    combined
    | collect
    | map{
        assert it.size() > 0, "In pairSamplesByStrategy with sample selection strategy ${strategy}, no samples were paired, likely due to 'same' or 'different' filter."
    }

    return combined
}