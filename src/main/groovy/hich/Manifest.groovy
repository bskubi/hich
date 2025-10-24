package hich
class Manifest {
    static Closure FastqAlign = { id, fq, aln, aln_idx_dir, aln_idx_pre, conf, cpus ->
            new FastqAlign(id, fq, aln, aln_idx_dir, aln_idx_pre, conf, cpus)
        }
}