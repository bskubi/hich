include {Chromsizes} from './chromsizes.nf'
include {GenomeReference} from './genomeReference.nf'
include {AlignerIndex} from './alignerIndex.nf'
include {FragmentIndex} from './fragmentIndex.nf'
include {emptyOnLastStep} from './extraops.nf'

/*
    Returns true if all conditions are true:
        - key is in the hashmap
        - value is truthy
        - value's string representation has 1+ characters
*/
def truthyString = {
    map, key ->
    def valueExists = map.get(key)
    def stringHasPositiveLength = map.get(key).toString().trim()

    valueExists && stringHasPositiveLength
}

/*
    Arguments:
        map:    the HashMap (i.e. a single sample HashMap)
        keyMap: a HashMap where keys are a list of attributes of map, values are lists of return values
        defaultVal: a default value to return
    Returns:
        The return value associated with a key in keyMap that's uniquely valid in map, or defaultVal if
        there are no keys in keyMap uniquely valid in map.

    'valid' here means that truthyString evaluates to true, which means that the key exists in the map,
    its value is truthy, and its string representation has a length >= 1

    Example:
        If keyMap contains 3 key:value pairs: [fastq1, fastq2]:fastq, [sambam]:sambam, [pairs]:pairs
        and has the defaultValue null, then it will return:
            'fastq' if map.fastq1 and map.fastq2 are valid and map.sambam and map.pairs are not valid
            'sambam' if map.sambam is valid and none of map.fastq1, map.fastq2, or map.pairs are valid
            'pairs' if map.pairs is valid and none of map.fastq1, map.fastq2, or map.sambam are valid
            null otherwise.
*/
def keySwitch = {
    map, keyMap, defaultVal ->

    // Get a set of the attributes of map that we're searching for, i.e. ['fastq1', 'fastq2', 'sambam', 'pairs']
    def keySet = keyMap.keySet().flatten().toSet()

    /* Iterate through the items in keyMap (i.e. keys = ['fastq1', 'fastq2'], val = 'fastq')
       Determine whether 'map' has all members of 'keys' (i.e. map.fastq1 and map.fastq2 exist),
       and does NOT have any of the members of the other keys in keyMap (i.e. there's no map.sambam or map.pairs)

       findResult extracts the first non-null value returned by the transformation
    */
    def result = keyMap.findResult {
        keys, val ->

        // The list of keys that must be present and must NOT be present in order to return 'val'
        def mustHave = keys                             // Every item that IS in 'keys'
        def mustNotHave = keySet.findAll{!(it in keys)} // Every item that is NOT in 'keys' but IS in 'keySet'
        
        // Check that the key:value pair is valid or NOT valid in map (i.e. in the sample HashMap) as appropriate
        if (mustHave.every{truthyString(map, it)}
            && mustNotHave.every{!truthyString(map, it)}) {
            return val
        }

        // Return null if the key:value pair isn't valid, which will be ignored by findResult
        return null
    }

    // result will be empty if no key:value pairs were uniquely valid.
    // If there was no uniquely present key in keySet, return the default value instead
    result ?: defaultVal
}

// Arguments: a single sample HashMap
// Returns: ["datatype":datatype] by either extracting 'datatype' if specified or inferring it.
// Infers sample.datatype based on whether the hashmap has uniquely valid values for fastq1+fastq2, sambam, or pairs
def getDatatype = {
    sample ->

    // We will either extract or infer the value of datatype.
    // Only fastq, sambam, and pairs datatypes are handled by Hich at present.
    def datatype = ""
    def handledDatatypes = ["fastq", "sambam", "pairs"]

    if (truthyString(sample, "datatype")) {
        // Use the given value of datatype if it is explicitly specified by the user
        datatype = sample["datatype"]
    }
    else {
        // If the user did not specify an explicit datatype, infer it from the
        // unique presence/absence of truthyString values for fastq1+fastq2, sambam, or pairs
        def inferredDatatype = keySwitch(sample,
                [
                    ["fastq1", "fastq2"]:"fastq",
                    ["sambam"]:"sambam",
                    ["pairs"]:"pairs"
                ],
                "unknown")
        datatype = inferredDatatype
    }
    
    // Check for an invalid datatype, either because the datatype could not be inferred
    // or because the user explicitly specified a datatype Hich do not handle.
    def invalidDatatype = (datatype == "unknown" || !(datatype in handledDatatypes))

    // Throw an exception if the datatype inference failed
    if (invalidDatatype) {
        throw new IllegalArgumentException(
            "During sample ingestion, datatype was ${datatype}. " +
            "Hich expects the datatype to be in ${handledDatatypes}. " +
            "Either specify this explicitly in the 'datatype' column of the sample file, " +
            "or ensure that every sample row has non-blank values for exactly one of the following sets of columns: " +
            "fastq1 and fastq2, sambam, or pairs. Problematic sample: ${sample.id}.\n\n${sample}"
            )
    }

    datatype
}

workflow Setup {
    take:
        samples
    
    main:
        
        
        samples
            | map {
                sample ->
                
                ////////////////////////////////////////////////////
                // Input validation for the sample.
                // 1. Set up a humid run if specified
                // 2. Add default params
                // 3. Set sample datatype
                // 4. Set sample id
                // 5. If the sample id matches any ConfigMap "ids" lists under
                //    params, update the sample
                // 6. Convert string paths to data files to file objects

                ///////////////////////////////
                // Humid run
                // Set the n_reads parameter to take only the first n_reads of the input file
                // Currently only works on fastq
                if (params.get("humid")) {
                    useDefaultHumidReadCount     = params.humid instanceof Boolean
                    defaultHumidReadCount        = params.general.humidDefault
                    sampleSpecificHumidReadCount = params.get("humid")

                    n_reads = useDefaultHumidReadCount ? defaultHumidReadCount : sampleSpecificHumidReadCount
                    if (n_reads) {
                        sample += ["n_reads": n_reads]
                    }
                }

                ///////////////////////////////////////////////////////
                // Default params
                // In params.defaults, we can set default values to be associated
                // with sample HashMap keys if the key is not specified by the user.
                // This checks each of the keys listed in params.defaults.
                // If the key is not in the sample HashMap (i.e. not specified by the user)
                // then set it to the default value.
                params.defaults.each {
                    key, defaultVal ->

                    key in sample ? null : (sample += [key:defaultVal])
                }

                ////////////////////////////////////
                // Ensure that sample has a 'datatype' entry (fastq, sambam, or pairs)
                // This refers to the original input datatype.
                sample += ["datatype": getDatatype(sample)]

                /*
                    Ensure that the condition, biorep, and techrep are either null
                    or that they have at least one non-whitespace character.
                */
                sample.condition ? sample += ["condition":sample.get(condition).toString().trim()] : null
                sample.biorep ? sample += ["biorep":sample.get(biorep).toString().trim()] : null
                sample.techrep ? sample += ["techrep":sample.get(techrep).toString().trim()] : null

                /////////////////////////////////
                // If an id is not explicitly given by the user for the sample,
                // Or is just whitespace, create one based on its condition, biorep, and techrep
                if (!truthyString(sample, "id")) {
                    // Potential identifiers to build the id string
                    // trim them to ensure the user doesn't accidentally leave whitespace in the cell
                    // and have this create weird problems down the line.
                    condition = sample.get('condition').trim() ? sample.condition : null
                    biorep = sample.get('biorep').trim() ? sample.biorep : null
                    techrep = sample.get('techrep').trim() ? sample.techrep : null

                    // Ensure that at least one identifier is present
                    if (!condition && !biorep && !techrep) {
                        throw new Exception(
                            "Samples with no explicit 'id' must be given at least a condition, biorep, or techrep " +
                            "as identifiers to construct an id value. The constructed id value will be the non-null"
                            "identifiers in order of condition, biorep, and techrep, separated by _, as in condition_biorep_techrep."
                        )
                    }

                    // Select non-null identifiers and separate by _
                    // if all are present, the id string is condition_biorep_techrep
                    identifiers = [condition, biorep, techrep].removeAll {[null]}
                    idString = identifiers.join("_").toString()

                    // Update the sample with the constructed id string.
                    sample += ["id":idString]
                }

                //////////////////////////////////////////////////////////////
                // ConfigMaps
                // It's possible to define sub-ConfigMaps within params that are
                // sample id-specific, which offsets another way to define settings
                // for particular sample subsets. This can be useful for defining
                // sample settings that are arrays or non-flat structures like HashMaps.
                // These will replace any default values or values in the sample file.
                params.each {
                    k, bundle ->

                    // Overwrite previous params with id-specific params
                    // Priority is to later-specified params in nextflow.config
                    bundleExists = bundle
                    bundleIsConfigMap = bundle.getClass() == nextflow.config.ConfigMap
                    bundleHasIDList = bundle.get("ids")
                    bundleParams = bundle.findAll{it.key != "ids"}
                    sampleID = sample.id.toString()
                    
                    if (bundleExists && bundleIsConfigMap && bundleHasIDList) {
                        bundleIDStrings = bundle.ids.collect{it.toString()}
                        bundleAppliesToSample = sampleID in bundleIDStrings

                        // Apply the update to samples whose id matches the "ids"
                        // list for the params bundle.
                        bundleAppliesToSample ? (sample += bundleParams) : null
                    }
                }
                
                ///////////////////////////////////////////////////////////////
                // Convert string paths to data files into file objects
                ingest = ["fastq1", "fastq2", "sambam", "pairs"]
                ingest.each {
                    ingest ->
                    // This should give an error if the file does not exist
                    // or if there is no data file specified for the sample.
                    ingest in sample ? sample += [ingest:file(sample[ingest])] : null
                }
                
                // Emit the updated sample
                sample
            }
            | GenomeReference   // Create any needed reference files for the sample
            | Chromsizes
            | AlignerIndex
            | FragmentIndex
            | set{samples}

    // Validate that all samples have distinct truthyString values of 'id'.
    // Collect them into a single list of HashMaps.
    // Extract the list of id values
    // Ensure that every id is a truthyString and that the size of the set
    // equals the size of the list to confirm they are all unique.
    samples
        | collect   // Collect all samples in the channel to a list
        | map {
            sampleHashMaps ->

            //sampleHashMaps is a list of all the sample HashMaps in the channel.

            // Extract all the raw ids from the sample HashMaps
            ids = sampleHashMaps.findAll{sample -> sample.id}

            // Get the set of unique truthyString ids from the extracted ids
            truthyStringIDs = ids.findAll{id -> truthyString(id)}

            // Ensure that there's an equal number of unique truthyString ids
            // as in the original id list, ensuring that all ids are unique
            // and are truthyStrings.
            allUniqueIDs = ids.size() == truthyStringIDs.toSet().size()
            
            if (!allUniqueIDs) {
                // If there are any non-truthyString ids or duplicates, throw an exception.
                idCounts = ids.countBy{it}
                throw Exception(
                    "Not all samples have unique ids or some ids are blank. The count of each id is: ${idCounts}.\n" +
                    
                    "Where id values come from:\n" +
                    "You can specify unique id values in the sample file in the 'id' column. " +
                    "Alternatively, if no 'id' is specified for a given sample, it will be constructed " +
                    "from the values of any present identifiers 'condition', 'techrep', and 'biorep', separated by underscores.\n" +
                    
                    "Diagnosing and fixing the problem:\n" +
                    
                    "1. Hich's initial default value of techrep is 1 and this will be applied to every condition and biorep where an " +
                    "explicit value of 'techrep' is not given, so if you have 2+ bioreps for a particular condition and did not specify " +
                    "their ids or techrep numbers, then this would create a conflict in setting the 'id'.\n" +
                    
                    "2. When Hich merges techreps to bioreps, it generates an id condition_biorep. When it merges bioreps to techreps, " +
                    "it generates an id equal to the shared condition name alone. If you removed defaults.techrep or defaults.biorep " +
                    "and did not specify an explicit techrep or biorep label for a particular sample, then a merge step could produce " +
                    "a new sample that conflicts with it. You can fix this by specifying an explicit id for the ingested condition or biorep sample.\n" +
                    
                    "3. Check the sample file to ensure you did not accidentally specify 2+ samples with the same id value.\n" +

                    "4. All ids must be strings with at least one non-whitespace character. Hich does substantial input validation to attempt" +
                    "to ensure this will always be true so it typically should not be a problem. If all else fails, specify unique per-sample ids " +
                    "and condition, biorep, and techrep values explicitly."
                    )
            }
        }

    samples = emptyOnLastStep("Setup", samples)

    emit:
        samples
}
