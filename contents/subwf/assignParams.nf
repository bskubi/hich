include {MakeMissingChromsizes} from './chromsizes.nf'
include {TryDownloadMissingReferences} from './genomeReferences.nf'
include {MakeMissingIndex} from './alignerIndex.nf'
include {MakeMissingDigest} from './makeDigest.nf'
include {emptyOnLastStep} from './extraops.nf'

def validKey = {
    map, key ->

    map.get(key) && map.get(key).toString().trim()
}

def keySwitch = {
    map, keyMap, defaultVal ->

    def keySet = keyMap.keySet().flatten().toSet()
    def result = keyMap.findResult {
        keys, val ->

        def truthy = keys
        def falsey = keySet.findAll{!(it in keys)}
        
        if (truthy.every{validKey(map, it)} && falsey.every{!validKey(map, it)}) {
            return val
        }
        return null
    }
    result ?: defaultVal
}

def getDatatype = {
    sample ->

    ["datatype":keySwitch(sample,
              [
                ["fastq1", "fastq2"]:"fastq",
                ["sambam"]:"sambam",
                ["pairs"]:"pairs"
              ],
              sample.get("datatype"))]
}

workflow AssignParams {
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
                if (params.get("humid")) {
                    def n_reads = params.humid instanceof Boolean ? params.general.humidDefault : params.humid
                    if (n_reads) {
                        sample += ["n_reads": n_reads]
                    }
                }

                ///////////////////////////////////////////////////////
                // Default params
                params.defaults.each {
                    k, v ->

                    
                    !(k in sample) ? sample += [(k):v] : null
                }

                ////////////////////////////////////
                // Datatype
                sample += getDatatype(sample)

                /////////////////////////////////
                // Id
                if (!validKey(sample, "id")) {
                    sample += ["id":"${sample.condition}_${sample.biorep}_${sample.techrep}".toString()]
                }

                //////////////////////////////////////////////////////////////
                // ConfigMaps
                params.each {
                    k, bundle ->

                    // Overwrite previous params with id-specific params
                    // Priority is to later-specified params in nextflow.config
                    if (bundle && bundle.getClass() == nextflow.config.ConfigMap && bundle.get("ids")) {
                        id = sample.id.toString()
                        bundle_ids = bundle.ids.collect{it.toString()}

                        if (bundle.containsKey("ids") && id in bundle_ids) {
                            update = bundle.clone()
                            update.remove("ids")
                            sample += update
                        }
                    }
                }
                
                ///////////////////////////////////////////////////////////////
                // File objects

                ["fastq1", "fastq2", "sambam", "pairs"].each {
                    // This should give an error if the file does not exist
                    // or if there is no data file specified for the sample.
                    it in sample ? sample[it] = file(sample[it]) : null
                }
                
                sample
            }
            | TryDownloadMissingReferences
            | MakeMissingChromsizes
            | MakeMissingIndex
            | MakeMissingDigest
            | set{samples}

    samples = emptyOnLastStep("setup", samples)

    emit:
        samples
}
