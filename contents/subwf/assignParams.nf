include {MakeMissingChromsizes} from './chromsizes.nf'
include {TryDownloadMissingReferences} from './genomeReferences.nf'
include {MakeMissingIndex} from './alignerIndex.nf'
include {MakeMissingDigest} from './makeDigest.nf'

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
                
                if (params.containsKey("humid") && !sample.get("n_reads")) {
                    sample += ["n_reads": params.general.humid.n_reads]
                }

                params.defaults.each {
                    k, v ->

                    // Add default params to sample hashmaps if not present
                    !(k in sample) ? sample += [(k):v] : null
                }

                sample += getDatatype(sample)
            } | map {
                sample ->
                
                if (!validKey(sample, "id")) {
                    sample += ["id":"${sample.condition}_${sample.biorep}_${sample.techrep}".toString()]
                }
                sample
            } | map{
                sample ->

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
                
                sample
            } | map {
                sample ->
                
                ["fastq1", "fastq2", "sambam", "pairs"].each {
                    // This should give an error if the file does not exist
                    // or if there is no sample in the file. Also we should
                    // be able to figure out the datatype from this.
                    it in sample ? sample[it] = file(sample[it]) : null
                }
                sample
            }
            | TryDownloadMissingReferences
            | MakeMissingChromsizes
            | MakeMissingIndex
            | MakeMissingDigest
            | set{samples}


    emit:
        samples
}