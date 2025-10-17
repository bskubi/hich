**Sample:** Data files and processing parameters represented as key-value pairs in a hashmap within the Hich Nextflow pipeline. As Hich progresses, each sample accumulates *sample attributes*, including paths pointing to results of incremental data *processing steps*.
**Sample attribute:** A single key:value pair in a sample.
**Input sample:** A *sample* that was input by the user (i.e. via a *sample file*)
**Constructed sample:** A *sample* constructed by Hich based on other *samples*, configuration, and *sample attributes* (i.e. as a result of an *aggregation plan* or *analysis plan*)
**Active sample:** The set of *input* and *constructed samples* that have not been filtered out by Hich.
**Active sample set:** The complete collection of *active samples*.
**Option:** A *sample attribute* that tells Hich how to configure one of more *processing steps*
**Power option:** A high-level *sample attribute* that Hich converts to one or more *direct options.*
**Direct option:** A low-level *sample attribute* that Hich directly passes to a specific CLI tool on a single *processing step*.
**Hardcoded option:** A low-priority default option that is overridden by any supplied *user options* that supply a value to the same *direct option.*
**User option:** A user-supplied option for a specific sample and processing step that override any *hardcoded options*.
**opts attribute:** A *sample attribute* in which the value is a hashmap containing *power options* and *direct options* used to configure a specific *processing step.* They are called *opts attributes* because they end in the suffix `_opts` (i.e. align_opts, matrix_opts, etc).
**Command option:** A hashmap within an *opts attribute* that contains direct options for a specific CLI command. Example: 
**Processing step:** Roughly corresponding to a set of Nextflow processes that all conceptually perform the same task (i.e. alignment, where we don't care about the specific aligner).
**Command:** String that invokes a tool to perform an individual processing step on a single sample.
**Aggregation level:** Condition, biological replicate, or technical replicate.
**Analysis plan:** A set of *sample selection strategies* plus a set of parameters to apply to every *sample* matching the combined *sample selection strategy.*
**Composite sample selection strategy:** A set of *sample selection strategies* used to identify a subset of the *active sample set*.
**Sample selection strategy:** A named set of *sample attribute* keys and values.
**Sample file:** A columnar file (i.e. tsv/csv) with a header in which rows define individual *input samples* and columns define *sample attributes*.
**Read:** A single fastq read (or read pair).
**Alignment:** A single alignment from a read (a read can produce more than one alignment).
**Pair:** A single contact between two positions parsed from alignments originating from one read.
***Cis* pair:** Intrachromosomal contact
***Trans* pair:** Interchromosomal contact

