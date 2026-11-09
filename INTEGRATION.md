# Independent R and edgeR integration check

Install the package with its test extra and provide Rscript with edgeR and limma. From an extracted source release:

```sh
python -m pip install '.[test]'
python integration/check_r_reference.py --output new_r_check
```

The script fails, rather than silently skips, if R or a required R package is absent. It uses only synthetic counts and a new output directory. It compares Python pseudobulk against independently coded R `rowsum`, and the descriptive natural-log1p CPM contrast against a separate R expression. A different synthetic fixture exercises the real edgeR adapter and an independent R edgeR implementation under the same design.

The two estimators are deliberately not equated: Python descriptive contrasts are not edgeR log2FC estimates. Seeds and numerical tolerances are fixed in the script before execution. Local execution on macOS ARM completed with R 4.5.2, edgeR 4.8.2 and limma 3.66.0: 12 donors/67 cells in the aggregation fixture; 200 genes in the edgeR fixture. Pseudobulk sums matched exactly; descriptive contrasts matched within 1e-11; edgeR logFC, logCPM, F, PValue and FDR matched within rtol 1e-10 / atol 1e-12.

These results establish numerical/interface agreement on synthetic inputs. They do not establish full SEA-AD reanalysis validity, sex-specific effects or causal biology. The script creates runtime diagnostics only in the user-selected directory; inspect those before sharing.
