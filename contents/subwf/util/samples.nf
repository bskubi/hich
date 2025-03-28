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