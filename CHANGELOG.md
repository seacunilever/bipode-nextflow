# seacunilever/bipode-nextflow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0 - [28/08/2026]

### `Added`

- **Multiple batch keys supported**
  - The `prepare-inputs` command now supports specifying multiple batch keys by repeating the `--batch-key` argument.
  - This enables more flexible control of batch correction and experimental grouping when multiple sources of technical variation are present.
  - Users should ensure that combinations of batch keys are linearly separable. Poorly separated batch structures may lead to reduced sampling efficiency during Bayesian model fitting.

- **Simplified binomial likelihood implementation**
  - The Bayesian concentration-response model has been updated to use a single binomial likelihood formulation across all count ranges.
  - Approximate likelihood formulations previously used for very high and very low count data have been removed.
  - Internal testing demonstrated that these approximations were unnecessary and that the simplified implementation provides equivalent inference while reducing model complexity and easing future maintenance.

### `Fixed`
- Update areas failing with strict syntax used in current Nextflow versions

## 1.0.0 - [25/03/2026]

### `Added`
- Change repo name to 'bipode-nextflow'
- Name change to bipode-httr Python package and version bump to 1.0.0

        load "nft-bam@0.6.0"
        load "nft-utils@0.0.9"

## 0.5.0 - [17/10/2025]

### `Added`
- Version bump for bifrost-httr Python package

### `Added`
- Version bump for bifrost-httr Python package

## 0.4.3 - []

### `Added`
- Compatibility with 0.4.3 bifrost-httr Python package

### `Fixed`
- Change references from 'seqera-services' to 'seacunilever'

### `Dependencies`

### `Deprecated`
- Remove unneeded action `bifrost/.github/workflows/awsfulltest.yml`

## 0.4.2 - [31/07/2025]

Initial release of seqera-services/bifrost, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
