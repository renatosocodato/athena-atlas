# athena-atlas

Donor-aware atlas analysis. Author: Renato Socodato. Version 0.1.0rc3, prepared for review; not published.

## Install and run

```sh
python -m pip install .
athena-atlas example --output example_output
athena-atlas run --input input.csv --output analysis_output
python -m pytest tests -q
```

Atlas additionally needs `--metadata metadata.csv`. Targetome accepts JSON instead of CSV. Output directories must be new; existing data are never overwritten by the public CLI. Analysis runs stage their outputs and materialize the analysis directory only after input validation and calculation succeed. Examples are deliberately synthetic and do not validate biological or clinical performance. See INPUTS.md for units, grouping and assumptions.

The core API is in `athena_atlas.core`. No network, cloud account or private directory is required. `tests/` checks scientific boundaries and invalid inputs as well as known reference answers.


## Interpretation

No therapeutic efficacy, causal disease mechanism, biomarker qualification or clinical performance follows from executing these tools. Synthetic examples establish executable behavior only. Full external-engine and full-cohort scientific revalidation remains separate from software packaging checks.

The descriptive API uses the sum of supplied genes as its count-per-million denominator and does not adjust for sex, batch or disease stage. See INPUTS.md before interpreting a filtered-gene matrix.

## Independent integration reference

See [INTEGRATION.md](INTEGRATION.md) for the executed R/edgeR synthetic reference and its limits.


## Scope of this release

The package distributes the documented public API and command-line workflows. The historical grant workflows, project data and unvalidated experimental stages are outside this release.
