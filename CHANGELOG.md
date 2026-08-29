# seacunilever/bipode-nextflow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0 - [29/08/2026]

### `Added`

- **Multiple batch keys supported**

  - The `prepare-inputs` command now supports specifying multiple batch keys by repeating the `--batch-key` argument.
  - This enables more flexible control of batch correction and experimental grouping when multiple sources of technical variation are present.
  - Users should ensure that combinations of batch keys are linearly separable. Poorly separated batch structures may lead to reduced sampling efficiency during Bayesian model fitting.

- **Simplified binomial likelihood implementation**
  - The Bayesian concentration-response model has been updated to use a single binomial likelihood formulation across all count ranges.
  - Approximate likelihood formulations previously used for very high and very low count data have been removed.
  - Internal testing demonstrated that these approximations were unnecessary and that the simplified implementation provides equivalent inference while reducing model complexity and easing future maintenance.

### `Changed`

- Renamed the remaining `BIFROST`/`bifrost` references to `BIPODE`/`bipode`, covering test fixtures, example output, documentation images and logo assets. Earlier `CHANGELOG` entries are left as-is, since they record the names in use at the time.

### `Fixed`

- Update areas failing with strict syntax used in current Nextflow versions
- Regenerated the `assets/test_data/minimal` fixtures against bipode-httr 1.1.0. The committed fixtures still used the pre-1.1.0 data schema (`n_treatment_batch`, flat `batch_index`, `low`/`high` count buckets), so the isolated module tests failed even though full pipeline runs — which generate their intermediates live — passed.
- Repaired the CI linting workflows: removed the orphaned `linting_comment` workflow, which waited on a `linting-logs` artifact that no longer exists since the `nf-core lint` job was stripped from `linting.yml`, and corrected the repository/branch guards in `fix_linting` and `awsfulltest`, none of which matched `seacunilever/bipode-nextflow`.
- Fixed the pre-commit violations (Prettier formatting, trailing whitespace, missing final newlines) that were failing the `nf-core linting` workflow.
- Corrected invalid defaults in `nextflow_schema.json`: the `seed` default contradicted the `null` default in `nextflow.config` and broke `nf-core pipelines lint`, and `batch_mode` was documented as `all` while the config sets `batch`.

### `Removed`

- Dropped stale test fixtures and assets that nothing referenced: the pre-1.1.0 `Nitrofurantoin_HepG2_*` files, the orphaned compiled model `assets/model/BIFROST_HTTr_beta_logistic_batch`, the duplicate `fix-linting.yml` workflow, and the accidentally committed `nextflow.config.amltmp`.

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
