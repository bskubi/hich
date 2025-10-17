# Hich Attributes
**Reliable:** Ensures the user gets the behavior they expect.
**Versatile:** Adapts to diverse use cases.
**Powerful:** Handles large, complex analysis tasks.
**Clear:** Easy to figure out how to get it working.
# Hich Components
+ **Interface** (config, params, profiles, sampleFile)
+ **Orchestration** (Nextflow, containers)
+ **Preprocessing** (alignment, filtering and contact matrix generation)
+ **Analysis** (feature calling and QC)
# Hich Capabilities
| **Component**     | **Reliable (Expected Behavior)**                                    | **Versatile (Diverse Cases)**                                | **Powerful (Handles Scale)**                                     | **Clear (Easy to Understand)**                                        |
| ----------------- | ------------------------------------------------------------------- | ------------------------------------------------------------ | ---------------------------------------------------------------- | --------------------------------------------------------------------- |
| **Interface**     | Provides input validation with informative errors.                  | Allows all internal tools to be fully configured.            | Enables precision control with analysis plans.                   | Offers a zero-install experience with declarative sample attributes.  |
| **Orchestration** | Uses a containerized, industry-standard workflow manager.           | Is compatible with a wide variety of computing environments. | Scales efficiently to handle massive datasets and sample counts. | Makes pipeline steps, inputs, and outputs predictable and observable. |
| **Preprocessing** | Is built on standard bioinformatics tools or well-tested new tools. | Accepts a wide variety of assays and data formats.           | Can merge and split replicates and cell types.                   | Auto-generates resource files to minimize user setup work.            |
| **Analysis**      | Employs validated, Hi-C-specific analysis methods.                  | Facilitates sensitivity analysis and data subsampling.       | Performs comprehensive Hi-C analysis and QC.                     | Generates interactive, detailed, GUI-based QC reports.                |
![tests](Tests/TEST_MATRIX.csv)

# Test productivity
#### Why not use [GitHub Actions](https://docs.github.com/en/actions/get-started/understand-github-actions) to run Hich tests?
* I looked into this on 2025-10-14. **GitHub actions is only free if using standard, free runners. The ones they offer only have 4 CPUs, 16GB RAM, and 14GB storage, which is too little for some of the tests we'll need to run.** [source](https://docs.github.com/en/actions/how-tos/write-workflows/choose-where-workflows-run/choose-the-runner-for-a-job#standard-github-hosted-runners-for-public-repositories)
* I'm not clear on whether it would be possible or helpful to use self-hosted runners on ARC, but this seems like a can of worms. I'd rather develop 1+ SLURM scripts that test Hich on ARC.
# nf-test issues
* nf-test shards by skipping tests, which could interact oddly with obsolete snapshot detection, but I don't know if this is the case.
* 5xing the shards for the fastest runs only decreased runtime by a factor of 2, possibly due to outlier slow tests.
* Is there a way to rerun only tests that failed on the last run?
* nf-test does not let you view the input variables to the function (AFAIK), so it's not possible to use it to test functions meant to modify variables passed by reference -- only functions that return values.
# Basics

[**nf-test**](https://www.nf-test.com/) is the automated test suite for Hich.