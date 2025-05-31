# seqera-services/bifrost

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# Bifrost

Bifrost is a Nextflow pipeline for analyzing high-throughput transcriptomics data (HTTr) to identify concentration-dependent responses and estimate points of departure (PoDs).

![Bifrost Pipeline Report](assets/images/pipeline_report.png)

## Introduction

The Bifrost model (Bayesian inference for region of signal threshold) is a statistical model for analysis of HTTr concentration-response data. The model is designed to infer a point-of-departure (PoD) from a concentration-response dataset. The PoD is an estimate of the minimum effect concentration of the test substance for the experimental conditions under which the data were produced. PoDs are estimated as probability distributions.

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

For more details about the output files and reports, please refer to the [output documentation](docs/output.md). An [example output report](docs/examples/BIFROST_input_Nitrofurantoin_HepaG2_full.html.zip) is available for Nitrofurantoin treatment in HepaG2 cells.

## Credits

Bifrost was originally written by [Joe Reynolds](https://github.com/JoeReynolds257) and [Mark Liddell](https://github.com/mark-liddell). [Jonathan Manning](https://github.com/pinin4fjords) later updated the workflow structure and migrated to nf-core standards.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
