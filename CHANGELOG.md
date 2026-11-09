# 0.1.0rc4

- Preserve literal CSV identifiers and reject duplicate headers and ambiguous metadata. Add a supported edgeR command with explicit independent donor metadata; exclude the historical H5AD/scaffold workflow.

# 0.1.0rc3

- Reject malformed/empty schemas and convert numeric input before aggregation.
- Stage CLI analysis output so invalid inputs leave no analysis directory.
- Refine optional dependencies to this package and separate heavy engines.
- 39 source-level tests passed on Python 3.14; final installed-distribution checks are recorded separately.
- Validate donor/cell/gene identifiers, integer count range and two-condition contrasts. Explain supplied-gene normalization and unadjusted biological scope.

# 0.1.0rc1

Initial private review candidate with installable API, explicit inputs and synthetic examples.

## 0.1.0rc3

Added an executable independent R aggregation reference and real edgeR adapter/reference comparison on synthetic counts. No private atlas data required.
