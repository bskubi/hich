package hich

class Hub {
    static Closure FastqAlign = { id, fastq1, fastq2, aligner_index_dir, config_fastq_align, cpus ->
        new FastqAlign(id, fastq1, fastq2, aligner_index_dir, config_fastq_align, cpus)
    }

    static Closure BamParsePairs = { id, bam, chromsizes, config_bam_parse_pairs, cpus, memory ->
        new BamParsePairs(id, bam, chromsizes, config_bam_parse_pairs, cpus, memory)
    }

    static Closure PairsLabel = { id, pairs, fragment_index, config_pairs_label ->
        new PairsLabel(id, pairs, fragment_index, config_pairs_label)
    }

    static Closure PairsSelect = { id, pairs, config_pairs_select, cpus ->
        new PairsSelect(id, pairs, config_pairs_select, cpus)
    }

    static Closure PairsMerge = { id, pairs, config_pairs_merge, cpus ->
        new PairsMerge(id, pairs, config_pairs_merge, cpus)
    }

    static Closure PairsDedup = { id, pairs, config_pairs_dedup, cpus ->
        new PairsDedup(id, pairs, config_pairs_dedup, cpus)
    }

    static Closure PairsCoolBin = { id, pairs, chromsizes, config_pairs_cool_bin ->
        new PairsCoolBin(id, pairs, chromsizes, config_pairs_cool_bin)
    }

    static Closure CoolCoarsenAddNorm = { id, cool, config_cool_coarsen_addnorm, cpus ->
        new CoolCoarsenAddNorm(id, cool, config_cool_coarsen_addnorm, cpus)
    }

    static Closure PairsHicBinCoarsenAddNorm = { id, pairs, chromsizes, config_hic_bin_coarsen_addnorm, cpus, memory ->
        new PairsHicBinCoarsenAddNorm(id, pairs, chromsizes, config_hic_bin_coarsen_addnorm, cpus, memory)
    }

    static Closure MatrixHiCRep = { id, matrix, command_hich_matrix_hicrep ->
        new MatrixHiCRep(id, matrix, command_hich_matrix_hicrep)
    }
}