# Pipeline Parameters

This document describes the parameters used in the Bifrost pipeline. Parameters are grouped into different sections based on their functionality.

For information about testing the pipeline and validating parameter inputs, see the [test guide](test_guide.md).

## Table of Contents

- [Input/Output Options](#inputoutput-options)
- [Bifrost-specific Options](#bifrost-specific-options)
- [Report Options](#report-options)
- [Generic Options](#generic-options)

## Input/Output Options

These parameters define where the pipeline should find input data and save output data.

| Parameter               | Type           | Description                                                                                                                                                                                                     | Required | Default |
| ----------------------- | -------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- | ------- |
| `input`                 | file path      | Path to input file(s). Choose ONE of these input modes:<br>1. A CSV file containing sample information (raw data input mode)<br>2. One or more pre-prepared JSON files in Bifrost format (prepared input mode). | Yes      | -       |
| `meta_mapper`           | file path      | Path to the metadata mapper YAML file that defines column mappings. Only required when processing raw data inputs.                                                                                              | No       | -       |
| `counts`                | file path      | Path to the counts CSV file containing probe counts. Only required when processing raw data inputs.                                                                                                             | No       | -       |
| `substances_cell_types` | file path      | Path to YAML file containing test substances and cell types to analyze. Only required when processing raw data inputs.                                                                                          | No       | -       |
| `batch_mode`            | string         | Way in which to collect probes before passing for analysis. Options:<br>- `all`: collect all probes in one tar file<br>- `batch`: group probe files by batch                                                    | No       | `batch` |
| `outdir`                | directory path | The output directory where the results will be saved. Use absolute paths for Cloud infrastructure.                                                                                                              | Yes      | -       |

## Bifrost-specific Options

Parameters specific to the Bifrost pipeline functionality.

| Parameter                  | Type    | Description                                                                                                                                                                                                        | Required | Default             |
| -------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------- | ------------------- |
| `seed`                     | integer | Random seed for reproducible results. If not set, results will not be reproducible between runs.                                                                                                                   | No       | -                   |
| `model_file`               | string  | Path to the Stan model file (without .stan extension)                                                                                                                                                              | No       | -                   |
| `n_cores`                  | integer | Number of CPU cores to use for processing                                                                                                                                                                          | No       | 4                   |
| `max_nodes`                | integer | Maximum number of nodes to use in Azure Batch (only applicable for Azure profile)                                                                                                                                  | No       | 1                   |
| `batch_key`                | string  | Field to use as batch key in the BIFROST model                                                                                                                                                                     | No       | "Exposure plate ID" |
| `min_percent_mapped_reads` | number  | Minimum percentage of mapped reads required (0-100)                                                                                                                                                                | No       | 50                  |
| `min_num_mapped_reads`     | integer | Minimum number of mapped reads required                                                                                                                                                                            | No       | 100000              |
| `min_avg_treatment_count`  | number  | Minimum average treatment count required                                                                                                                                                                           | No       | 5                   |
| `test_probes`              | integer | Number of probes to sample for testing (optional)                                                                                                                                                                  | No       | -                   |
| `precompile_model`         | boolean | Whether to pre-compile the Nextflow script before execution. Set to 'true' on shared file systems to compile the model once and share it across instances. Set to 'false' when no shared file system is available. | No       | -                   |

## Report Options

Parameters for customizing the MultiQC report generation.

| Parameter                     | Type    | Description                                                                          | Required | Default    |
| ----------------------------- | ------- | ------------------------------------------------------------------------------------ | -------- | ---------- |
| `report_timepoint`            | string  | Exposure duration within experiment                                                  | No       | "24 hours" |
| `report_conc_units`           | string  | Concentration units for test substance (`uM`, `ugml-1`, or `mgml-1`)                 | No       | "uM"       |
| `report_interactive_plots`    | boolean | Force interactive plots (may be faster for large datasets)                           | No       | true       |
| `report_n_fold_change_probes` | integer | Number of most up/down regulated probes to show                                      | No       | 2          |
| `report_cds_threshold`        | number  | Concentration-Dependency Score threshold for filtering probes (0-1)                  | No       | 0.5        |
| `report_n_lowest_means`       | integer | Number of lowest mean PoD probes to show                                             | No       | 10         |
| `report_n_pod_stats`          | integer | Number of probes to include in PoD statistics table                                  | No       | 100        |
| `report_plot_height`          | integer | Height of concentration-response plots in pixels                                     | No       | 400        |
| `report_pod_vs_fc_height`     | integer | Height of PoD vs Fold Change plot in pixels                                          | No       | 600        |
| `report_no_cds_threshold`     | boolean | Do not filter probes by CDS threshold in summary tables and lowest mean PoDs section | No       | false      |

## Generic Options

Less common options for the pipeline, typically set in a config file.

| Parameter          | Type    | Description                                                                                                                     | Required | Default |
| ------------------ | ------- | ------------------------------------------------------------------------------------------------------------------------------- | -------- | ------- |
| `publish_dir_mode` | string  | Method used to save pipeline results to output directory. Options: `symlink`, `rellink`, `link`, `copy`, `copyNoFollow`, `move` | No       | "copy"  |
| `monochrome_logs`  | boolean | Do not use coloured log outputs                                                                                                 | No       | false   |
| `validate_params`  | boolean | Validate parameters against the schema at runtime                                                                               | No       | true    |

### Hidden Options

These options are typically not modified by users but are available for advanced use cases:

- `version`: Display version and exit
- `pipelines_testdata_base_path`: Base URL for pipeline test dataset files
- `trace_report_suffix`: Suffix for trace report filename (default: date and time)
