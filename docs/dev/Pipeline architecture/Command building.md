Each *processing step* has a `buildCmd` function in its `functions.nf` file. The `buildCmd` method receives *sample attributes* and returns a *command*. It starts with *hardcoded options*, overrides them with any conflicting *user options*, and maps *direct options* with multiple different keys (i.e. `-o` and `--output`) to the same final *direct option* key.

Create tests for each `buildCmd` function to ensure it's written correctly, which will also enable me to refactor it with confidence #TestPlan

Refactor `buildCmd` so that it outputs a list of values combined in the process, which would simplify testing it. #Refactoring
