include {coalesce} from './reshape.nf'

def label(map, lbl) {
    // Return whether map.lbl contains a non-empty string, used below to determine
    // if the techrep, biorep and condition keys are present and specified
    return map?[(lbl)]?.toString()?.length() > 0
}

/*
    Determine if the techrep, biorep, condition fields are uniquely specified

    For isSingleCell, the user has to specify a cellBarcodeField (which tag in a sam/bam
    file holds the cell barcode in order to extract it to the pairs cellID field)
    or has to specify isSingleCell (which can be used for pairs files where the cellID field
    is already extracted).
*/
def isTechrep(map) {return label(map, "techrep") && label(map, "biorep") && label(map, "condition")}
def isBiorep(map) {return (!label(map, "techrep")) && label(map, "biorep") && label(map, "condition")}
def isCondition(map) {return (!label(map, "techrep")) && (!label(map, "biorep")) && label(map, "condition")}
def isSingleCell(map) {return map.cellBarcodeField || map.isSingleCell}

/*
    Hich depends on each sample having a unique sample.id attribute for joining process results to the
    appropriate sample hashmap. A legible name is also convenient for troubleshooting. If the user
    wants to let Hich build unique ids automatically, they should specify unique conditions, bioreps and techreps
    and not use the _ character in order to ensure that all ids will be unique. The aggregateProfileName is also
    included because new copies of the input samples are produced for each aggregateProfile.
*/
def makeID(map, columns) {
    
    if (columns) {
        map = coalesce(map, true)
    }

    return (
        map.subMap("condition", "biorep", "cell", "techrep", "aggregateProfileName").values().join('_')
    )
}

def aggregateLevelLabel(sample) {
    // Return a magic string for the aggregateLevel
    if (isTechrep(sample)) return "techrep"
    if (isBiorep(sample)) return "biorep"
    if (isCondition(sample)) return "condition"
    return "unknown"
}

