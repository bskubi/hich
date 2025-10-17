# `_opts maps` in sample attributes

Config options for individual tasks is stored in sample attributes suffixed with `_opts`. These `_opts` attributes are hashmaps. There are two types: **direct options**, passed unmodified to the underlying CLI tool, and **power options**, higher-level sample config that Hich converts to direct options. 

In SelectPairs

Power options:
+ `select_pairs_opts.filters`: optional
	+ `.discardSingleFrag` -- defaults to true if fragmentIndex is present, false otherwise
	+ `.keepPairTypes` -- defaults to `["UU", "UR", "RU"]`
	+ `.onlyCis`
	+ `.onlyTrans`
	+ `.minDistFR`
	+ `.minDistRF`
	+ `.minDistFF`
	+ `.minDistRR`
	+ `.custom`
	+ `.chroms`

