[GitHub issue](https://github.com/bskubi/hich/issues/38)
```
SmokeTest_MergeBiorepsToConditions  
SmokeTest_DifferentialLoops  
SmokeTest_Select_buildCmdEmpty  
projectDirSampleFile  
singleEntrySampleFile
```

Test: projectDirSampleFile and singleEntrySampleFile
+ Test that an empty input channel with sample_file_single_entry.tsv passes
+ Problem is that the sample file has been removed/eliminated entirely, and it's not just a rename -- there is no sample_file_single_entry.tsv in the test assets.
+ Solution: figure out what this is supposed to test and recreate an appropriate sample file.
	+ Created a sample file that matches assert conditions in `assets/sampleFile/single_entry.tsv`
+ This test now passes.
+ singleEntrySampleFile was essentially redundant, so I removed it.


Test: SmokeTest_Select_buildCmdEmpty
+ Tests buildCmd function in SelectPairs/functions.nf
+ There are now 5 params for buildCmd, so both buildCmd and buildCmd_Empty should fail
+ The Select buildCmd method now takes a fragment index as input and uses this to determine whether to discard fragments unless the user specifies otherwise in select_pairs_opts
+ There still appear to be mismatched snake_case and mixedCase variable names, such as fragmentIndex and select_pairs_opts in buildCmd for SelectPairs.
+ Solution: Test with or without a fake fragment index being passed in, as well as the override capability.
+ SmokeTest_Select_buildCmdEmpty expected that there were no defaults, but this is no longer true. I removed it as we probably need to rethink our overall test strategy rather than do ad-hoc test fixes for this issue.

Test: SmokeTest_DifferentialLoops
+ ERROR ~ Cannot enable more than one container engine -- Choose either one of: docker, apptainer
+ Uses -c docker.config, which enables docker. The problem is that nextflow.config in the root project directory enables apptainer.
+ Running again gives a new error: `.command.sh: line 2: diff_mustache: command not found`
+ This is because diff_mustache is not installed, but needs to be called via python -m. I made this change, and the test now passes.

Test: SmokeTest_MergeBiorepsToConditions
+ Error is that it now returns 14 samples instead of the expected 13.
+ There are 10 input sample files, with 2x bioreps and 4x conditions, merged as follows:
	+ 1.1+1.2
	+ 2.1+2.1
+ 5 conditions observed, expected 4, correct number of techreps and bioreps
	+ We have two condition=4, both with same id.
	+ The idea is that we already have a replicate with condition = 4 at aggregateLevel condition, while the other condition = 4 samples should all have merging disabled.
	+ Removing the condition=4, id=4, aggregateLevel="condition" sample causes the correct number of conditions to be produced.
	+ Bad: 4_1_skipbiorepmerge is being merged to a condition
	+ Good: 4_1_skipmerge is not being merged to a condition
	+ Bad: 4_1_nomerge is not being merged to condition
	+ Good: 4_1_1 is not being merged to condition
	+ First, looks like 4_1_nomerge is being labeled with mergeBiorepToCondition: false, when it should be more general. Currently, it's redundant with 4_1_skipbiorepmerge.
	+ Looks like mergeBiorepToCondition: false is not being respected by the pipeline.
	+ We renamed the skip steps as:
		+ skipTechrepMerge
		+ skipBiorepMerge
		+ skipMerge
	+ Test now passes

```
                    [id: "4_1_1"] + [condition: "4", biorep: "1", techrep: "1", aggregateLevel: "techrep", pairs: pairs, latestPairs: pairs],
                    [id: "4"] + c4 + [aggregateLevel: "condition", pairs: pairs, latestPairs: pairs],
                    [id: "4_1_nomerge"] + b1 + c4 + bioreps + [mergeBiorepToCondition: false],
                    [id: "4_1_skipmerge"] + b1 + c4 + bioreps + [skipMerge: true],
                    [id: "4_1_skipbiorepmerge"] + b1 + c4 + bioreps + [mergeBiorepToCondition: false]
```
# Final Status:
+ Passes all --tag 1to3s tests
# Changes
+ Process DIFFERENTIAL_LOOPS now calls `python -m diff_mustache` instead of `diff_mustache`
+ Removed `SmokeTest_Select_buildCmdEmpty`
* Added `assets/sampleFile/single_entry.tsv` and used this for test `projectDirSampleFile`
* Check techrep, biorep, condition sizes before total size in SmokeTest_MergeBiorepsToConditions 
* Updated `skip.*Merge` attributes in SmokeTest_MergeBiorepsToConditions

# Future directions
+ Process LOOPS also needs to be called with `python -m mustache` instead of just `mustache`
+ Test the slower nf-tests
