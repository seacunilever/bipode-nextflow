# Bifrost Nextflow workflow

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

Bifrost is a bioinformatics pipeline for processing and analyzing concentration response data. The workflow is designed to be portable across different execution environments (local, HPC, cloud providers) and takes in pre-formatted json data files for one cell line/chemical (referred to as one dataset).

## Usage

> [!NOTE]
> If you are new to Nextflow, please refer to [this page](https://www.nextflow.io/docs/latest/getstarted.html) on how to set-up Nextflow. Make sure to [test your setup](https://www.nextflow.io/docs/latest/getstarted.html#testing-the-installation) with `-profile test` before running the workflow on actual data.

First, prepare your input data in the required format. Then, you can run the pipeline using:

```bash
nextflow run bifrost \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --counts counts.csv \
   --outdir <OUTDIR>
```

For more details and further functionality, please refer to the [usage documentation](docs/usage.md) and the [output documentation](docs/output.md).

## Pipeline output

For more details about the output files and reports, please refer to the [output documentation](docs/output.md).

## Credits

Bifrost was originally written by [Joe Reynolds](https://github.com/JoeReynolds257) and [Mark Liddell](https://github.com/mark-liddell). [Jonathan Manning](https://github.com/pinin4fjords) later updated the workflow structure and migrated to nf-core standards.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
