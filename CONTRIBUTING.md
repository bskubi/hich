# Structure

Hub

# YAML config

- Channel items must be sortable for Nextflow to ensure reproducibility and
to enable caching. This requires that objects have homogeneous structure.
In particular, keys like 'id' should be all strings or all integers. May be
best to make all keys strings.

# Naming conventions

- Follow nf-core naming conventions

# How to avoid Nextflow's pitfalls

- Variables declared using `def` in a process script/stub block shadow variables declared in input/output and prevents them from taking on values.
- Variables declared in Nextflow map operator closures must be declared using def. If not, cryptic concurrency bugs result if the same variable name is reused in a later map closure.
- Emitting classes seems to result in strange concurrency bugs that I haven't solved.

# Hich toolkit notes

## Alignment

- The `-p` flag should be used exclusively for interleaved FASTQ files, as only single end alignments and corrupt pairs result from using it on paired-end reads.