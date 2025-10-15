[github issue](https://github.com/bskubi/hich/issues/40)

SmokeTest_LabelMatrixPlans
+ doesn't find id2_plan2, aggregateLevel=biorep
+ This may be due to changes in how aggregation plans are specified. The test puts the values in matrices.plan1 or matrices.plan2
+ Currently, nf-test is not outputting the actual contents when using assertContainsInAnyOrder.
+ the assertContainsAnyOrder appears to be complaining that there are unexpected samples in the output, not expected samples missing from the output.
+ The issue was that the output samples lacked a "resolutions" parameter. Unclear on if that's an issue or not.
+ This is likely because rather than putting the options into the sample directly, they are stored under an attribute matrix_opts (which includes the planName and sampleSelectionStrategy).

# Changes
+ Expect that matrix plans store params under sample matrix_opts attribute in SmokeTest_LabelMatrixPlans