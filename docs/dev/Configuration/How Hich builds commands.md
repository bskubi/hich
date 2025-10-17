The CLI command for each process is built in a *processing step*-specific `buildCmd` method in `Step/functions.nf`. This step starts with low-priority *hardcoded options* that it overrides with *user options*, some of which are obtained by converting *power options* to *direct options*.  


*Opts attributes* may contain a mix of *power options* and *direct options*. *Power options* are converted into *direct options* by Hich for the individual sample while the CLI command is being built. Opts attributes may be composed in part of nested *power opts attributes* (i.e. align_opts may contain bwameth_opts or bwa_mem_opts)
### align_opts
**Power options:**
+ `bwa_mem_opts` (hashmap): *Direct options* for `bwa mem` and `bwa-mem2 mem`.
+ `bwameth_opts` (hashmap): *Direct options* for `bwameth`.
+ `filterMAPQ` (boolean): If true, then the *sample attribute* `minMapq` is supplied to the aligner as a quality threshold to report alignments.
+ Note: `align_opts` itself does not contain *direct options* for aligners, those must be supplied via `align_opts.bwa_mem_opts` or `align_opts.bwameth_opts` depending on the value of the `aligner` *sample attribute*.
### ingest_pairs_opts
**Power options:**
**Direct options:**
### parse_to_pairs_opts
**Power options:**
+ `.pairtools_parse2_opts`
+ `.hich_pairs_sql`
+ `.pairtools_sort`
**Direct options:**

### select_pairs_opts
**Power options:**
+ `.filters (hashmap)`: *Power options* for `pairtools select`
	+ `.discardSingleFrag` (`boolean`): Discard pairs if both ends map to the same restriction fragment (default: `true` if sample has `fragmentIndex` , `false` otherwise).
	+ `.keepPairTypes` (`list[str]`): List of pair types to retain (default: `["UU", "RU", "UR"]`).
	+ `.onlyCis` (`boolean`): If true, keep only *cis pairs*.
	+ `.onlyTrans` (`boolean`): If true, keep only *trans pairs*.
	+ `.minDistFR` (`integer`): Minimum end-to-end distance to keep on FR-orientation pairs.
	+ `.minDistRF`: Minimum end-to-end distance to keep on RF-orientation pairs.
	+ `.minDistFF`: Minimum end-to-end distance to keep on FF-orientation pairs.
	+ `.minDistRR`: Minimum end-to-end distance to keep on RR-orientation pairs.
	+ `.custom`: Arbitrary string passed to *pairtools select*.
	+ `.chroms`: List of chromosome names that are written to a `.bed` file and provided to `--chrom-subset`.
**Direct options:**
+ Any key:value pairs in `select_pairs_opts` other than `.filters` are passed as *direct options* to `pairtools select`

### tag_restriction_fragments_opts
**Power options:**
+ `hich_pairs_map_ends_opts (hashmap)`: Contains direct options passed to `hich pairs map-ends`
**Direct options:**
### matrix_opts
**Power options:**
+ `resolutions (list[integer])`:
+ `juicer_tools_pre (hashmap)`:
+ `cooler_zoomify_opts (hashmap)`:
+ `cooler_cload_pairs (hashmap)`:
+ `cooler_balance_opts (hashmap)`:
**Direct options:**