# Input contract

Counts CSV: unique nonmissing cell index; unique nonempty genes as columns; finite, nonnegative integer counts. Numeric strings are converted to integer counts before aggregation. Empty matrices, invalid counts and values beyond the supported exact integer range are rejected. Metadata CSV: exactly the same cell index in the same order, with donor, condition and sex columns. Donor labels must identify independent biological donors. Each donor must have one condition and sex annotation. Empty annotations are rejected.

`athena-atlas run` requires `--metadata metadata.csv`, and optional `--treatment` / `--control` select the two conditions. Subset to those two conditions explicitly; additional conditions are rejected. At least two donors per condition are required, but this is an execution minimum and does not establish adequate biological precision.

The public descriptive contrast normalizes by the sum of supplied genes per donor, then computes natural-log1p counts per million and a difference of donor means. Supply the same full gene/count universe; a selected gene panel changes the denominator and therefore the estimand. The value is not a log2 fold change. Bootstrap intervals resample independent donors within condition. Two donors cannot support a well-resolved interval simply by increasing cell count or bootstrap draws.

The descriptive API does not adjust for sex, batch or Braak stage, fit paired/repeated-donor designs, or establish causal disease effects. CSV cell and donor identifiers are read literally; duplicate headers and surrounding whitespace are rejected rather than inventing genes or splitting donor identities.

## Separate edgeR contrast

`athena-atlas edger --counts donor_counts.csv --metadata donor_metadata.csv --route female --output edge_result` uses a gene-by-sample raw count matrix (gene IDs in the first column). The metadata must contain one row per sample and independent donor, with `sample_id`, `donor_id`, `sex_route`, and `diagnosis`. Count columns must match metadata sample IDs exactly. In the selected route, diagnosis is `control` or `ad`, with at least two independent donors per group as an execution minimum. Optional `braak_numeric` must be finite; a varying covariate enters the design, and non-estimable designs fail explicitly. Install R/edgeR/limma separately. The output is a quasi-likelihood contrast for that design, with log2 fold changes and statistical descriptors; it does not establish a causal disease effect or validate any particular cohort. The command does not read or round normalized H5AD data.
