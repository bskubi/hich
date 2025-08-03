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
            attribute, acceptableValues ->
            sample.get(attribute) in acceptableValues
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
    /* Get all samples having matching values of strategy.same
    */
    return samples
        | map{tuple(it.subMap(strategy.get("same", [])), it)}
        | groupTuple
        | map{it[1]}
}

def pairSamplesByStrategy(samples, strategy) {
    /*
        The diffloops workflow has to define a way to pair up samples, especially
        by enforcing that certain attributes are the same or different.

        samples -- channel of sample hashmaps
        strategy -- a composite sample selection strategy [attribute: permittedValues]
            two possible attributes are
                same: a list of attributes which must be the same to pair two samples
                different: a list of attributes which must differ to pair two samples
            other attributes are used for filtering individual samples
    */

    // Ensure there's no conflict between "same" and "different"
    def same = strategy.get("same", [])
    def different = strategy.get("different", [])
    def sameAndDifferent = same.intersect(different)
    if (!sameAndDifferent.isEmpty()) {
        System.err.println("Warning: In filterSamplesByStrategy, comparisons on ${sameAndDifferent} are required to be same and different, so no result is obtained")
        return channel.empty()
    }

    // Filter individual samples before forming pairs of samples
    def filtered = filterSamplesByStrategy(samples, strategy)

    // Obtain pairs of samples matching the "same" and "different" criteria
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