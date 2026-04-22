# Hich: A Modular Pipeline for Hi-C and DNA Methylation Across Unimodal, Co-assay, Bulk, and Single-Cell Formats

Hich is a unified processing pipeline for bulk and single-cell Hi-C and DNA methylation assays. Where existing tools bottleneck on single-purpose designs, hardcoded parameters, and restricted query capabilities, Hich provides a modular architecture:
+ **Unification:** Processes both unimodal and multi-omic coassays within a single framework.
+ **Extensibility:** Replaces hardcoded parameters with Python-based User-Defined Functions (UDFs) for dynamic filtering, reshaping, and quality control.
+ **Scalability:** Leverages TileDB for out-of-core condition, cell, and spatial slicing of base-pair resolution contacts and per-read cytosine states.
+ **Interoperability:** Reshapes multidimensional array slices for immediate integration with downstream aggregation and analysis tools.