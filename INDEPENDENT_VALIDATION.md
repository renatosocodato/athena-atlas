# Independent implementation and scientific checks

Independent R aggregation exactly reproduced donor pseudobulk counts for 12 synthetic donors and 67 cells. The descriptive natural-log1p CPM contrast agreed within 1e-11. The packaged edgeR wrapper and a separate R implementation agreed for logFC, logCPM, F, PValue and FDR in 200 synthetic genes at rtol 1e-10 / atol 1e-12. The rc4 public command reproduced the reference with macOS R 4.5.2 / edgeR 4.8.2, preserving literal sample IDs. The earlier rc3 adapter was also compared on Ubuntu ARM64.

This establishes implementation equivalence on the tested counts, not disease effects or biological validity of a complete human atlas. The descriptive API is a different estimator from edgeR. See integration/ for the portable reference check.

These are internal computational checks against independent implementations, public inputs and real engines. No external peer review is claimed. Full execution details are in VALIDATION.md.

## Revision after independent candidate review

Preserve literal CSV identifiers and reject duplicate headers and ambiguous metadata. Add a supported edgeR command with explicit independent donor metadata; exclude the historical H5AD/scaffold workflow. Hosted execution records for the exact distributions are tracked separately.
