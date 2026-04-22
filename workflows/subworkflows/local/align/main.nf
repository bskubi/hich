workflow Align {
    take:
    samples
    is3DGenome
    ALIGN_PROC
    

    main:
    samples
        | filter{ it.fastqs }
        | map{ tuple(it.ID, it.fastqs, it.indexDir, it.indexPfx, is3DGenome, it.isInterleaved) }
        | ALIGN_PROC
        | set{result}

    emit:
    result
}