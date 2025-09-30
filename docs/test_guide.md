# Testing

This document describes how to test the Bifrost pipeline, including both manual testing for your setup and automated testing with nf-test.

## Manual Testing

Before running the pipeline with your own data, it's recommended to test your setup using the provided test profile. The test profile uses a small dataset and minimal computational resources to verify that the pipeline runs correctly in your environment.

### Quick test

To run a quick test with Docker:

```bash
nextflow run seacunilever/bifrost -profile test,docker --outdir test_results
```

To run a quick test with Singularity:

```bash
nextflow run seacunilever/bifrost -profile test,singularity --outdir test_results
```

To run a quick test with Conda:

```bash
nextflow run seacunilever/bifrost -profile test,conda --outdir test_results
```

Alternatively, you can install the dependencies directly on your host system:

```bash
pip install bifrost-httr
# You'll also need to install cmdstan - see https://pypi.org/project/bifrost-httr/
nextflow run seacunilever/bifrost -profile test --outdir test_results
```

### Test profiles available

The pipeline provides a test profile:

- `test`: Minimal test using a small subset of data (5 probes, ~10 minutes runtime)

### What the test does

The test profile:

- Uses minimal test data from `assets/test_data/minimal/`
- Processes a small dataset with 5 probes across multiple concentrations
- Tests both Paracetamol and Nitrofurantoin treatments in HepG2 cells
- Runs with reduced computational resources (2 cores, 1 node)
- Generates a test report to verify all pipeline steps work correctly

### Expected test output

If the test runs successfully, you should see:

- Pipeline completion message
- Generated HTML report in the `test_results/report/` directory
- Compressed results in `test_results/compressed_results/`
- Pipeline execution reports in `test_results/pipeline_info/`

The test should complete in approximately 10 minutes depending on your system.

## Automated Testing

# seacunilever/bifrost: Testing

## Introduction

This document describes the testing framework and procedures for the Bifrost pipeline. The pipeline uses [nf-test](https://nf-co.re/docs/nf-test/overview) for testing individual modules and the complete pipeline workflow, along with JSON Schema validation for input validation.

## Testing Framework Overview

The testing framework consists of several components:

1. **nf-test**: Primary testing framework for Nextflow workflows and modules
2. **Pipeline Schema**: Parameter validation using Nextflow schema
3. **CI/CD integration**: Automated testing via GitHub Actions
4. **Test sharding**: Parallel test execution for improved performance
5. **Linting**: Code quality checks using nf-core standards

## nf-test Configuration

### Global Configuration

The pipeline uses `nf-test.config` to define global test settings:

```groovy
config {
    // location for all nf-test tests
    testsDir "."

    // nf-test directory including temporary files for each test
    workDir System.getenv("NFT_WORKDIR") ?: ".nf-test"

    // location of an optional nextflow.config file specific for executing tests
    configFile "tests/nextflow.config"

    // ignore tests coming from the nf-core/modules repo
    ignore 'modules/nf-core/**/tests/*', 'subworkflows/nf-core/**/tests/*'

    // run all test with defined profile(s) from the main nextflow.config
    profile "test"

    // list of filenames or patterns that should trigger a full test run
    triggers 'nextflow.config', 'nf-test.config', 'conf/test.config', 'tests/nextflow.config', 'tests/.nftignore'

    // load the necessary plugins
    plugins {
        load "nft-bam@0.4.0"
        load "nft-utils@0.0.3"
    }
}
```

### Test Types

The pipeline includes two main types of tests:

#### 1. Module-Level Tests

Individual modules have their own test suites located in `modules/local/*/tests/main.nf.test`. These tests verify that each module:

- Executes successfully with valid inputs
- Produces expected outputs
- Handles edge cases appropriately

**Example module test structure:**

```groovy
nextflow_process {
    name "Test Process PREPARE_INPUTS"
    script "../main.nf"
    process "PREPARE_INPUTS"

    test("Should run without failures") {
        when {
            process {
                """
                input[0] = file("${projectDir}/assets/test_data/minimal/Example_Meta_Data.csv", checkIfExists: true)
                input[1] = file("${projectDir}/assets/test_data/minimal/sers_meta_data_mapper.yml", checkIfExists: true)
                input[2] = file("${projectDir}/assets/test_data/minimal/Example_Counts_5probes.csv", checkIfExists: true)
                input[3] = file("${projectDir}/assets/test_data/minimal/substances_cell_types.yml", checkIfExists: true)
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out.prepared_inputs).match()
        }
    }
}
```

#### 2. Pipeline-Level Tests

The global test suite is located in `tests/default.nf.test` and tests the complete pipeline workflow:

```groovy
nextflow_pipeline {
    name "Test pipeline with default settings"
    script "../main.nf"

    test("Params: prepared inputs") {
        when {
            params {
                outdir = "$outputDir"
                input = "${projectDir}/assets/test_data/minimal/BIFROST_input_Nitrofurantoin_HepG2.json"
                meta_mapper = null
                counts = null
                substances_cell_types = null
            }
        }

        then {
            def stable_name = getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
            assertAll(
                { assert workflow.success},
                { assert snapshot(workflow.trace.succeeded().size(), stable_name).match() }
            )
        }
    }
}
```

### Available Test Modules

The pipeline includes tests for the following modules:

- `COMPILE_STAN_MODEL`: Tests Stan model compilation
- `COMPRESS_OUTPUT`: Tests output compression functionality
- `CONC_RESPONSE_ANALYSIS`: Tests concentration-response analysis
- `CREATE_MULTIQC_REPORT`: Tests report generation
- `PREPARE_INPUTS`: Tests input data preparation
- `SPLIT_DATA`: Tests data splitting functionality

## Schema Validation

### Pipeline Schema

The main pipeline schema is defined in `nextflow_schema.json` and validates:

- Parameter types and ranges
- Required vs optional parameters
- Default values
- Parameter descriptions and help text

## CI/CD Integration

### GitHub Actions Workflows

The pipeline uses several GitHub Actions workflows for automated testing:

#### 1. nf-test Workflow (`.github/workflows/nf-test.yml`)

**Triggers:**

- Pull requests (excluding documentation changes)
- Releases
- Manual dispatch

**Features:**

- **Test sharding**: Automatically distributes tests across multiple CI jobs
- **Matrix testing**: Tests multiple Nextflow versions and profiles
- **Concurrency control**: Cancels outdated runs to save resources
- **Artifact collection**: Stores test results and logs

**Test Matrix:**

- **Profiles**: `docker`
- **Nextflow versions**: `24.04.2`, `latest-everything`
- **Sharding**: Up to 7 parallel shards

#### 2. Linting Workflow (`.github/workflows/linting.yml`)

**Components:**

- **Pre-commit checks**: Code formatting and basic linting
- **nf-core linting**: Pipeline-specific linting using nf-core standards
- **Release validation**: Additional checks for release branches

### Test Sharding

The pipeline uses intelligent test sharding to improve CI performance:

1. **Shard calculation**: The `get-shards` action analyzes changed files and calculates optimal shard distribution
2. **Parallel execution**: Tests are distributed across multiple CI jobs
3. **Load balancing**: Ensures even distribution of test execution time
4. **Resource optimization**: Reduces total CI runtime and resource usage

## Running Tests Locally

### Prerequisites

1. Install nf-test:

```bash
# Using conda
conda install -c bioconda nf-test

# Or download directly
curl -fsSL https://code.askimed.com/install/nf-test | bash
```

2. Install Nextflow (version 24.04.2 or later)
3. Install Docker, Singularity, or Conda (depending on your preferred execution method)

### Running All Tests

```bash
# Run all tests with Docker
nf-test test --profile docker

# Run all tests with Singularity
nf-test test --profile singularity

# Run all tests with Conda
nf-test test --profile conda
```

### Running Specific Tests

```bash
# Run pipeline-level tests only
nf-test test tests/

# Run specific module tests
nf-test test modules/local/prepare_inputs/tests/

# Run tests with specific tags
nf-test test --tag preprocessing
```

### Running Tests in CI Mode

```bash
# Run tests as they would run in CI
nf-test test --profile docker --ci --changed-since HEAD^ --verbose --tap=test.tap
```

### Test Output

nf-test produces several types of output:

1. **TAP format**: Machine-readable test results
2. **Snapshots**: Stored expected outputs for comparison
3. **Test logs**: Detailed execution logs
4. **Nextflow logs**: Pipeline execution details

## Test Data

### Test Data Location

Test data is stored in `assets/test_data/` with two main datasets:

- `minimal/`: Small dataset for quick testing (5 probes, ~562 samples)
- `full/`: Comprehensive dataset for thorough testing

### Test Data Structure

```
assets/test_data/
├── minimal/
│   ├── BIFROST_input_Nitrofurantoin_HepG2.json     # Pre-prepared JSON input
│   ├── Example_Meta_Data.csv                        # Sample metadata
│   ├── Example_Counts_5probes.csv                   # Expression counts
│   ├── sers_meta_data_mapper.yml                    # Meta data mapping
│   └── substances_cell_types.yml                    # Analysis configuration
└── full/
    ├── Example_Counts.csv                           # Full expression data
    ├── Example_Meta_Data.csv                        # Full metadata
    └── sers_meta_data_mapper.yml                    # Meta data mapping
```

### Snapshot Management

When test outputs intentionally change due to pipeline updates, you'll need to update the snapshots:

```bash
# Update snapshots after intentional changes
nf-test test --profile docker --updateSnapshot

# You may need to re-specify the test profile if using another
nf-test test --profile test,local --updateSnapshot
```

## Troubleshooting

### Common Issues

1. **Test failures after code changes**: Update snapshots with `--update-snapshot`
2. **Schema validation errors**: Check pipeline parameters against `nextflow_schema.json`
3. **CI timeout issues**: Consider test sharding or reducing test data size
4. **Module test failures**: Verify module inputs and outputs match test expectations

### Verbose Mode

```bash
# Run with verbose logging to see underlying Nextflow commands
nf-test test --profile docker --verbose
```

### Test Cleanup

```bash
# Clean test work directories
nf-test clean

# Remove old test outputs
rm -rf .nf-test/
```

## Integration with Development Workflow

### Pre-commit Testing

Before committing changes:

1. Run affected tests locally
2. Update snapshots if outputs have intentionally changed
3. Verify CI checks will pass

### Pull Request Testing

When creating pull requests:

1. All tests must pass in CI
2. New functionality should include appropriate tests
3. Test coverage should be maintained or improved
4. Documentation should be updated for new test procedures

This comprehensive testing framework ensures the reliability and maintainability of the Bifrost pipeline across different environments and use cases.
