# seqera-services/bifrost: Usage

## Introduction

This document describes how to use the Bifrost pipeline. The pipeline is designed to be portable across different execution environments (local, HPC, cloud providers) and takes in pre-formatted json data files for one cell line/chemical (referred to as one dataset).

## Prerequisites

- Linux (required for Nextflow, can be WSL2 https://learn.microsoft.com/en-us/windows/wsl/install)
- Nextflow version 21.04.0 or later (https://www.nextflow.io/docs/latest/getstarted.html)

## Input Modes

The pipeline supports two distinct input modes:

1. Raw data input - for processing new data, requires a samplesheet, counts file, and configuration files
2. Pre-prepared JSON input - for using previously prepared data, accepts pre-formatted JSON files directly

Choose the appropriate mode based on your needs:

- Use raw data input when starting with new experimental data
- Use pre-prepared JSON input when working with data that has already been formatted for Bifrost

### Pre-prepared JSON Input

To use pre-prepared JSON files:

```bash
nextflow run bifrost \
   -profile <docker/singularity/.../institute> \
   --input prepared_data.json \
   --outdir <OUTDIR>
```

The JSON files must follow the Bifrost format. Here's an example (truncated for readability):

```json
{
    "test_substance": "Nitrofurantoin",
    "cell_type": "HepG2",
    "probes": ["ACBD3_59", "CEP89_11010", "MPDU1_11661", "OAS3_90233", "TMEM183A_34069"],
    "counts": [
        [79, 141, 153, 249, 64, 159, 39, ...],  // Counts for probe 1
        [75, 240, 140, 109, 61, 150, 21, ...],  // Counts for probe 2
        [263, 310, 319, 438, 202, 255, 35, ...], // Counts for probe 3
        [17, 104, 78, 69, 18, 54, 15, ...],     // Counts for probe 4
        [117, 187, 174, 220, 135, 159, 41, ...] // Counts for probe 5
    ],
    "batch_index": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, ...],
    "concentration": [0.0192, 0.096, 0.48, 2.4, 12.0, 60.0, 300.0, ...],
    "n_treatment_batch": 5
}
```

The JSON format consists of:

- `test_substance`: Name of the substance being tested
- `cell_type`: Type of cells used in the experiment
- `probes`: Array of probe identifiers
- `counts`: 2D array where each inner array contains counts for one probe across all samples
- `batch_index`: Array of batch indices (1-based) for each sample
- `concentration`: Array of concentrations for each sample
- `n_treatment_batch`: Number of treatment batches in the dataset

When using pre-prepared JSON input:

- The meta_mapper, counts, and substances_cell_types parameters are not required
- You can provide multiple JSON files using the --input parameter

## Raw Data Input

### Input File Requirements

For raw data input, you need three files with specific structures and formats:

1. **Metadata CSV file** (`--input`) - Sample information and experimental metadata
2. **Counts CSV file** (`--counts`) - Gene expression count data matrix  
3. **Configuration YAML file** (`--substances-cell-types`) - Analysis configuration

### Metadata CSV File Structure

The metadata file must be a comma-separated CSV file with specific required columns. The pipeline validates this file and applies quality filters before analysis.

#### Required Columns

| Column                 | Type    | Description                                                       | Validation Rules                    |
| ---------------------- | ------- | ----------------------------------------------------------------- | ----------------------------------- |
| `SAMPLE_ID`            | string  | Unique identifier for each sample                                 | Must not contain spaces            |
| `CELL_TYPE`            | string  | The type of cell used in the experiment                          | Must not contain spaces            |
| `TEST_SUBSTANCE`       | string  | The substance being tested                                        | Must not contain spaces            |
| `CONCENTRATION`        | numeric | The concentration of the test substance                           | Must be non-negative (≥ 0)         |
| `NUM_MAPPED_READS`     | integer | Number of mapped reads                                            | Must be non-negative (≥ 0)         |
| `PERCENT_MAPPED_READS` | numeric | Percentage of mapped reads                                        | Must be between 0 and 100          |

#### Optional Columns

| Column                | Type    | Description                                               | Usage                                    |
| --------------------- | ------- | --------------------------------------------------------- | ---------------------------------------- |
| `TREATMENT_VESSEL_ID` | string  | ID of the treatment vessel                                | Used as batch key by default             |

Additional columns present in your metadata file are preserved but not actively used by the pipeline analysis.

#### Quality Filtering Criteria

The pipeline automatically filters samples based on the following criteria (configurable via parameters):

- **Minimum percentage of mapped reads**: Default 50% (`--min_percent_mapped_reads`)
- **Minimum number of mapped reads**: Default 100,000 (`--min_num_mapped_reads`)  
- **Minimum average treatment count**: Default 5.0 (`--min_avg_treatment_count`)

Samples not meeting these criteria are excluded from analysis.

#### Sample ID Requirements

- Sample IDs must be unique across the entire dataset
- Sample IDs cannot contain spaces
- Sample IDs should be consistent with column headers in the counts file

#### Example Metadata Structure

Based on the test data, here's an example of the metadata file structure:

```csv
SAMPLE_ID,CELL_TYPE,TEST_SUBSTANCE,CONCENTRATION,NUM_MAPPED_READS,PERCENT_MAPPED_READS,TREATMENT_VESSEL_ID
S_O5180393_HG2_NFUR_1,HepG2,Nitrofurantoin,0.0192,2857440,86.0,A18039301
S_M5180393_HG2_NFUR_2,HepG2,Nitrofurantoin,0.096,5710831,95.35,A18039301
S_K5180393_HG2_NFUR_3,HepG2,Nitrofurantoin,0.48,4481281,84.35,A18039301
S_I5180393_HG2_NFUR_4,HepG2,Nitrofurantoin,2.4,5654424,95.05,A18039301
S_G5180393_HG2_NFUR_5,HepG2,Nitrofurantoin,12.0,3290920,78.26,A18039301
S_E5180393_HG2_NFUR_6,HepG2,Nitrofurantoin,60.0,6389756,95.9,A18039301
S_C5180393_HG2_NFUR_7,HepG2,Nitrofurantoin,300.0,1538838,76.1,A18039301
S_B10180393_HG2_DMSO_0,HepG2,DMSO,0.0,4842380,95.87,A18039301
```

This example shows:
- **Dose-response series**: Nitrofurantoin at concentrations from 0.0192 to 300 μM
- **Control samples**: DMSO controls with concentration 0.0
- **Quality metrics**: Mapped reads ranging from ~1.5M to ~6.4M with mapping percentages 76-96%
- **Batch information**: All samples from the same treatment vessel (A18039301)

### Counts CSV File Structure

The counts file contains the gene expression count data in a matrix format where:

- **Rows represent probes/genes** 
- **Columns represent samples**
- **First column contains probe identifiers**
- **Remaining columns contain count values for each sample**

#### File Format Requirements

- Must be a comma-separated CSV file
- First column must contain unique probe identifiers  
- Column headers must match `SAMPLE_ID` values from metadata file
- Count values must be non-negative integers
- Missing values are not allowed

#### Structure and Validation

The pipeline validates that:
- The file has at least 2 columns (probe ID + at least one sample)
- All count values are non-negative  
- Sample IDs in column headers match those in the metadata file
- No missing or invalid count data

#### Example Counts File Structure

Based on the test data, here's the expected format:

```csv
Unnamed: 0,S_O5180393_HG2_NFUR_1,S_M5180393_HG2_NFUR_2,S_K5180393_HG2_NFUR_3,S_I5180393_HG2_NFUR_4,...
ACBD3_59,79,141,153,249,...
CEP89_11010,75,240,140,109,...
MPDU1_11661,263,310,319,438,...
OAS3_90233,17,104,78,69,...
TMEM183A_34069,117,187,174,220,...
```

**Key characteristics:**
- **Probe identifiers**: First column contains unique probe names (e.g., `ACBD3_59`, `CEP89_11010`)
- **Sample columns**: Each subsequent column represents one sample with its count data
- **Count values**: Integer expression counts for each probe in each sample
- **Matrix dimensions**: In the test data, 5 probes × 562 samples

#### Data Processing Notes

- The pipeline will log-transform concentrations for modeling
- Counts are split into "high" (>100) and "low" (≤100) categories for statistical modeling
- Batch effects are modeled using the specified batch key column

### Samplesheet input

For raw data input, you will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with a header row, containing the required columns as defined in the schema.

```bash
--input '[path to samplesheet file]'
```

### Example Samplesheet

Here's an example of a minimal samplesheet for testing Nitrofurantoin on HepG2 cells:

```csv
SAMPLE_ID,CELL_TYPE,TEST_SUBSTANCE,CONCENTRATION,NUM_MAPPED_READS,PERCENT_MAPPED_READS,TREATMENT_VESSEL_ID
S_O5180393_HG2_NFUR_1,HepG2,Nitrofurantoin,0.0192,2857440,86.0,A18039301
S_M5180393_HG2_NFUR_2,HepG2,Nitrofurantoin,0.096,5710831,95.35,A18039301
S_K5180393_HG2_NFUR_3,HepG2,Nitrofurantoin,0.48,4481281,84.35,A18039301
S_I5180393_HG2_NFUR_4,HepG2,Nitrofurantoin,2.4,5654424,95.05,A18039301
S_G5180393_HG2_NFUR_5,HepG2,Nitrofurantoin,12.0,3290920,78.26,A18039301
S_E5180393_HG2_NFUR_6,HepG2,Nitrofurantoin,60.0,6389756,95.9,A18039301
S_C5180393_HG2_NFUR_7,HepG2,Nitrofurantoin,300.0,1538838,76.1,A18039301
```

### Substances and Cell Types Configuration

You also need to provide a YAML file specifying which test substances and cell types to analyze. This file should be provided using the `--substances-cell-types` parameter:

```bash
--substances-cell-types '[path to substances_cell_types.yml]'
```

#### Configuration File Structure

The YAML file must contain the following sections:

```yaml
# Test substances to analyze
Test substances:
  - Nitrofurantoin
  - Paracetamol

# Cell types to analyze  
Cell types:
  - HepG2

Additional divider: N/A

Specific filters: null
```

#### Configuration Options

- **`Test substances`**: List of substances to analyze. Must match values in the `TEST_SUBSTANCE` column of your metadata file.

- **`Cell types`**: List of cell types to analyze. Must match values in the `CELL_TYPE` column of your metadata file.

- **`Additional divider`**: Optional field to further subdivide the analysis. Options:
  - Set to a column name from your samplesheet to create separate analyses for each unique value in that column
  - For example, if set to `TREATMENT_VESSEL_ID`, the pipeline will create separate analyses for each treatment vessel
  - Set to `N/A` to disable additional subdivision

- **`Specific filters`**: Optional dictionary to exclude specific values from analysis. Example:
  ```yaml
  Specific filters:
    TREATMENT_VESSEL_ID:
      - A18039301  # Exclude this treatment vessel
    CELL_TYPE:
      - HepG2      # Exclude this cell type
  ```

#### Example Configuration

Based on the test data, here's a complete example:

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

This configuration tells the pipeline to:
1. Analyze only Nitrofurantoin experiments
2. Include only HepG2 cell data  
3. Not use additional subdivision
4. Apply no specific exclusion filters

### Batch Key Configuration

The pipeline uses a batch key to group samples for statistical analysis. By default, it uses the `TREATMENT_VESSEL_ID` column, but you can change this using the `--batch-key` parameter:

```bash
--batch-key 'YOUR_COLUMN_NAME'
```

The batch key should be a column in your samplesheet that identifies groups of samples that were processed together (e.g., same plate, same experiment, etc.). This is used to account for batch effects in the statistical model.

#### Choosing an Appropriate Batch Key

Common batch key options include:
- `TREATMENT_VESSEL_ID`: Groups samples by treatment vessel/plate (default)
- `CELL_BATCH_ID`: Groups samples by cell culture batch
- `SEQUENCING_PLATE_ID`: Groups samples by sequencing plate
- `MEASUREMENT_DATE`: Groups samples by measurement date

Choose a batch key that represents the most relevant source of technical variation in your experiment.

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

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

### Batching Controls

The pipeline provides two modes for handling probe data batching through the `--batch-mode` parameter:

- `batch` (default): Groups probe files for analysis in smaller batches. This produces more intermediate files but uses less disk space overall when not using a shared file system.
- `all` (default): Collects all probes into a single tar file that is sent to all analysis processes. This reduces the number of file operations but requires more disk space when not using a shared file system.

Use `batch` unless you have good reason to believe file operations are a limiting factor in your infrastructure.

### Model Pre-compilation (shared file systems only)

For shared file systems, the pipeline offers a pre-compilation option for the Stan model through the `--precompile-model` parameter, allowing the model to be compiled once and used by all fitting processes. This has modest performance gains, but importantly prevents concurrent writes to the compiled model path in shared file systems (see below).

Recommended settings:

- Use `--precompile-model` when running on shared file systems (e.g., HPC clusters, Cloud with shared filesystem like Fusion)

  - Required: Must have a shared file system accessible to all processes
  - Benefits:
    - Single compilation shared across all processes
    - Prevents race conditions during compilation
    - Efficient use of shared storage
  - Note: Will fail without a shared file system due to CmdStanPy limitations

- Do not use `--precompile-model` when running on non-shared systems (e.g., cloud executors)
  - Each process compiles its own copy of the model
  - Works with process-specific local disks

Note: The pipeline will still work on shared file systems without this flag (using `stageInMode=copy`), but this negates the benefits of shared storage by forcing local copies of all inputs.

### Resuming

Add `-resume` when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

### Cleanup

Nextflow keeps all logs and files generated for runs in the work directory unless they are removed, so the workflow can be resumed.

## Report Generation Parameters

The pipeline generates an interactive MultiQC report that can be customized using various parameters. These parameters control the content, appearance, and behavior of the report.

### Basic Report Settings

| Parameter                    | Description                                                     | Default    |
| ---------------------------- | --------------------------------------------------------------- | ---------- |
| `--report_timepoint`         | Exposure duration in the experiment                             | "24 hours" |
| `--report_conc_units`        | Units for concentration values                                  | "uM"       |
| `--report_interactive_plots` | Enable interactive plot mode (may be faster for large datasets) | false      |

### Probe Selection and Filtering

| Parameter                       | Description                                                       | Default |
| ------------------------------- | ----------------------------------------------------------------- | ------- |
| `--report_cds_threshold`        | Minimum Concentration-Dependency Score (CDS) for filtering probes | 0.5     |
| `--report_n_fold_change_probes` | Number of most up/down regulated probes to show                   | 5       |
| `--report_n_lowest_means`       | Number of lowest mean PoD probes to show                          | 10      |
| `--report_n_pod_stats`          | Number of probes to include in PoD statistics table               | 100     |

### Plot Configuration

| Parameter                         | Description                                      | Default |
| --------------------------------- | ------------------------------------------------ | ------- |
| `--report_plot_height`            | Height of concentration-response plots in pixels | 400     |
| `--report_pod_vs_fc_height`       | Height of PoD vs Fold Change plot in pixels      | 600     |
| `--report_control_line_tolerance` | Tolerance for filtering similar control lines    | 0.02    |
| `--report_min_control_lines`      | Minimum number of control lines to show          | 2       |

### Performance Settings

| Parameter                             | Description                              | Default |
| ------------------------------------- | ---------------------------------------- | ------- |
| `--report_timeout`                    | Timeout in seconds for report generation | 300     |
| `--report_plots_force_flat_numseries` | Maximum number of series for flat plots  | 10000   |

### Example Usage

To generate a report with custom settings:

```bash
nextflow run seqera-services/bifrost \
    --input samplesheet.csv \
    --counts counts.csv \
    --substances_cell_types config.yml \
    --report_timepoint "48 hours" \
    --report_conc_units "ugml-1" \
    --report_interactive_plots \
    --report_cds_threshold 0.6 \
    --report_n_fold_change_probes 10 \
    --report_n_lowest_means 15 \
    --report_plot_height 500 \
    --report_pod_vs_fc_height 800
```

> [!NOTE]
> For large datasets, consider using `--report_interactive_plots` as it may significantly improve report generation performance. However, this will require JavaScript to be enabled in the browser to view the plots.

## Running on Cloud Platforms

Nextflow supports execution on various cloud platforms including AWS Batch, Azure Batch, and Google Cloud Batch. nf-core provides pre-configured profiles to simplify cloud setup, allowing for cost-effective processing using managed compute resources.

### Cloud Configuration Profiles

nf-core provides pre-configured profiles for popular cloud platforms:

- **AWS Batch**: [nf-core/configs/awsbatch](https://nf-co.re/configs/awsbatch/)
- **Azure Batch**: [nf-core/configs/azurebatch](https://nf-co.re/configs/azurebatch/)
- **Google Cloud Batch**: Available through standard Nextflow configuration

#### Organization-Specific Profile

This pipeline also includes a `unilever_azure` profile configured for Unilever's Azure infrastructure. This profile:

- Uses Unilever's Azure Container Registry and Batch account
- Configures auto-scaling pools with `Standard_F4s_v2` VMs
- Uses low-priority nodes for cost optimization (up to `max_nodes` parameter)
- Automatically creates and deletes compute pools
- Sets retry strategy (3 attempts) for handling preemptible node interruptions

To use this profile (requires Unilever Azure credentials):

```bash
nextflow run seqera-services/bifrost -profile unilever_azure --input samplesheet.csv --outdir results
```

### Getting Started with Cloud Execution

To run the pipeline on a cloud platform:

1. **Choose your cloud platform** and review the nf-core configuration profile
2. **Set up the required cloud resources** (compute environments, storage, etc.)
3. **Configure authentication** for your cloud provider
4. **Run the pipeline** with the appropriate profile:

```bash
# AWS Batch example
nextflow run seqera-services/bifrost -profile awsbatch --input samplesheet.csv --outdir results

# Azure Batch example  
nextflow run seqera-services/bifrost -profile azurebatch --input samplesheet.csv --outdir results
```

### Cloud Platform Documentation

For detailed setup instructions and configuration options, refer to the official documentation:

#### nf-core Cloud Configs
- [nf-core/configs](https://nf-co.re/configs/) - Browse all available institutional and cloud configurations
- [AWS Batch config](https://nf-co.re/configs/awsbatch/) - AWS-specific configuration and setup
- [Azure Batch config](https://nf-co.re/configs/azurebatch/) - Azure-specific configuration and setup

#### Nextflow Cloud Documentation
- [AWS Batch](https://www.nextflow.io/docs/latest/aws.html#aws-batch) - Official Nextflow AWS documentation
- [Azure Batch](https://www.nextflow.io/docs/latest/azure.html#azure-batch) - Official Nextflow Azure documentation  
- [Google Cloud Batch](https://www.nextflow.io/docs/latest/google.html#cloud-batch) - Official Nextflow Google Cloud documentation

### Cloud Execution Tips

- **Storage**: Use cloud-native storage solutions (S3, Azure Blob, Google Cloud Storage) for input and output data
- **Costs**: Consider using spot/preemptible instances for cost savings
- **Networking**: Ensure proper VPC/network configuration for security and performance
- **Monitoring**: Use cloud-native monitoring tools or Nextflow Tower for pipeline execution tracking

## Developer Notes

Many execution platforms are inefficient if a workflow tries to execute many short running processes. It can take more time to schedule and request resources for each small instance than bundling the short processes into a larger process task. Nextflow channel operators collate groups together inputs into batches which can run for longer with the short tasks themselves parallelised inside the process script which is what is done in `conc_response_analysis.py`.

The tens of thousands of probe pkl files from `split_data.py` were taking almost 2hrs to be transferred to the cloud blob storage work directory due to the large number of files. The solution implemented here compresses the pkl files into a `.tar.gz`, then each of the concentration response modelling processes then extract the files they need from it.

In order to facilitate inspection of dataset results as early as possible during running, `groupTuple()` is used with a `groupKey()` which is produced from the number of fit tar.gz files or each dataset. Doing this means Nextflow knows when all the of fits are available to run the compress process for the dataset rather than waiting until all the fits for all the datasets are complete.

## Useful Links

- Nextflow docs for Azure https://www.nextflow.io/docs/latest/azure.html
- Setup for Azure
  - https://shaunchuah.github.io/posts/setting-up-azure-with-nextflow
  - https://seqera.io/blog/nextflow-and-azure-batch-part-1-of-2/
  - https://techcommunity.microsoft.com/t5/healthcare-and-life-sciences/covid-variant-analysis-on-azure-using-nextflow-part-2/ba-p/3075741
- Batch pool scaling https://learn.microsoft.com/en-us/azure/batch/batch-automatic-scaling
- Pool scaling formula was taken from here https://github.com/Azure/doAzureParallel/blob/master/R/autoscale.R

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

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
