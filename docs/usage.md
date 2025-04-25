# seqera-services/bifrost: Usage

## Introduction

This document describes how to use the Bifrost pipeline. The pipeline is designed to be portable across different execution environments (local, HPC, cloud providers) and takes in pre-formatted json data files for one cell line/chemical (referred to as one dataset).

## Prerequisites

- Linux (required for Nextflow, can be WSL2 https://learn.microsoft.com/en-us/windows/wsl/install)
- Nextflow version 21.04.0 or later (https://www.nextflow.io/docs/latest/getstarted.html)

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with a header row, containing the required columns as defined in the schema.

```bash
--input '[path to samplesheet file]'
```

### Required Columns

The samplesheet must contain the following required columns:

| Column | Description |
|--------|-------------|
| `SAMPLE_ID` | Unique identifier for each sample. Must not contain spaces. |
| `CELL_TYPE` | The type of cell used in the experiment. Must not contain spaces. |
| `TEST_SUBSTANCE` | The substance being tested. Must not contain spaces. |
| `CONCENTRATION` | The concentration of the test substance (numeric) |
| `NUM_MAPPED_READS` | Number of mapped reads (numeric) |
| `PERCENT_MAPPED_READS` | Percentage of mapped reads (numeric) |

### Optional Columns

The following columns are optional but may be required depending on your analysis:

| Column | Description |
|--------|-------------|
| `TREATMENT_VESSEL_ID` | ID of the treatment vessel (used as batch key by default) |
| `EXPOSURE_TIME` | Duration of exposure (numeric) |

### Additional Requirements

- The pipeline will filter samples based on the following criteria:
  - Minimum percentage of mapped reads (default: 50%)
  - Minimum number of mapped reads (default: 100,000)
  - Minimum average treatment count (default: 5.0)
- Sample IDs must not contain spaces
- Numeric values should be provided as numbers, not strings

### Example Samplesheet

Here's an example of a minimal samplesheet for testing Nitrofurantoin on HepG2 cells:

```csv
SAMPLE_ID,CELL_TYPE,TEST_SUBSTANCE,CONCENTRATION,NUM_MAPPED_READS,PERCENT_MAPPED_READS,TREATMENT_VESSEL_ID,EXPOSURE_TIME
S_O5180393_HG2_NFUR_1,HepG2,Nitrofurantoin,0.0192,2857440,86.0,A18039301,24.0
S_M5180393_HG2_NFUR_2,HepG2,Nitrofurantoin,0.096,5710831,95.35,A18039301,24.0
S_K5180393_HG2_NFUR_3,HepG2,Nitrofurantoin,0.48,4481281,84.35,A18039301,24.0
S_I5180393_HG2_NFUR_4,HepG2,Nitrofurantoin,2.4,5654424,95.05,A18039301,24.0
S_G5180393_HG2_NFUR_5,HepG2,Nitrofurantoin,12.0,3290920,78.26,A18039301,24.0
S_E5180393_HG2_NFUR_6,HepG2,Nitrofurantoin,60.0,6389756,95.9,A18039301,24.0
S_C5180393_HG2_NFUR_7,HepG2,Nitrofurantoin,300.0,1538838,76.1,A18039301,24.0
```

### Substances and Cell Types Configuration

You also need to provide a YAML file specifying which test substances and cell types to analyze. This file should be provided using the `--substances-cell-types` parameter:

```bash
--substances-cell-types '[path to substances_cell_types.yml]'
```

Example `substances_cell_types.yml`:
```yaml
# Test substances to analyze
Test substances:
  - Nitrofurantoin

# Cell types to analyze
Cell types:
  - HepG2

Additional divider: N/A

Specific filters: null
```

This configuration tells the pipeline to analyze Nitrofurantoin on HepG2 cells. The following fields are available:

- `Test substances`: List of substances to analyze
- `Cell types`: List of cell types to analyze
- `Additional divider`: Optional field to further subdivide the analysis. If set to a column name from your samplesheet, the pipeline will create separate analyses for each unique value in that column. For example, if set to `TREATMENT_VESSEL_ID`, it will create separate analyses for each treatment vessel. Set to `N/A` to disable.
- `Specific filters`: Optional dictionary of filters to exclude specific values. For example:
  ```yaml
  Specific filters:
    TREATMENT_VESSEL_ID:
      - A18039301  # Exclude this treatment vessel
    CELL_TYPE:
      - HepG2      # Exclude this cell type
  ```

### Batch Key Configuration

The pipeline uses a batch key to group samples for statistical analysis. By default, it uses the `TREATMENT_VESSEL_ID` column, but you can change this using the `--batch-key` parameter:

```bash
--batch-key 'YOUR_COLUMN_NAME'
```

The batch key should be a column in your samplesheet that identifies groups of samples that were processed together (e.g., same plate, same experiment, etc.). This is used to account for batch effects in the statistical model.

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run seqera-services/bifrost --input ./samplesheet.csv --counts ./counts.csv --substances_cell_types ./substances_cell_types.yml --outdir ./results  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
