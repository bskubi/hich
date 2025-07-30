include {isExistingFile} from '../../util/files.nf'

workflow PrepIndex {
    take:
        samples
    main:
        samples
            | filter {it.datatype == "fastq"}
            | filter {!isExistingFile(it.alignerIndexDir)}
            | map{it.alignerIndexPrefix = it.alignerIndexPrefix ?: it.assembly; it}
            | map{tuple(it.genomeReference, it.genomeReference, it.alignerIndexPrefix)}
            | unique
            | set{result}
    emit:
        result
}