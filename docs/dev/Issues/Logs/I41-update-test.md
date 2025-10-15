[github issue](https://github.com/bskubi/hich/issues/41)

SmokeTest_SelectPairs_buildCmd
+ doesn't find expected `--output` component of command
+ it is now named `-o`, likely because of the multi-option homogenization I added 
+ test should just be updated to use the appropriate output

# Final status
nf-test test --tag SmokeTest_SelectPairs_buildCmd passes

# Changes
+ Test now expects -o rather than --output