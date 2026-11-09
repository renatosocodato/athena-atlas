"""Raw-count aggregation and descriptive contrasts at the independent donor level."""
import numpy as np
import pandas as pd


def pseudobulk(counts, metadata):
    """Sum integer counts per donor after validating the biological hierarchy.

    Numeric strings are accepted and converted before summation. Missing donor
    annotations, duplicate cell/gene identifiers, and inconsistent donor metadata
    are rejected. Output counts are signed 64-bit integers.
    """
    if not isinstance(counts, pd.DataFrame) or not isinstance(metadata, pd.DataFrame):
        raise ValueError('Counts and metadata must be data frames')
    if counts.empty or not counts.columns.is_unique or counts.columns.hasnans:
        raise ValueError('Nonempty counts and unique nonmissing gene identifiers required')
    if any(not str(gene).strip() for gene in counts.columns):
        raise ValueError('Empty gene identifier')
    if not counts.index.equals(metadata.index) or not counts.index.is_unique or counts.index.hasnans:
        raise ValueError('Unique, aligned, nonmissing cell indices required')
    if any(not str(cell).strip() or str(cell) != str(cell).strip() for cell in counts.index):
        raise ValueError('Cell identifiers must be nonempty and without surrounding whitespace')
    required = ['donor', 'condition', 'sex']
    if not set(required).issubset(metadata.columns) or not metadata.columns.is_unique:
        raise ValueError('Unique donor, condition and sex metadata columns required')
    if metadata[required].isna().any().any():
        raise ValueError('Missing donor metadata')
    for column in required:
        labels = metadata[column].astype(str)
        if labels.str.strip().eq('').any() or labels.ne(labels.str.strip()).any():
            raise ValueError('Donor metadata must be nonempty and without surrounding whitespace')
        if metadata[column].nunique() != labels.nunique():
            raise ValueError(f'Ambiguous mixed-type {column} identifiers')
    try:
        numeric = counts.apply(pd.to_numeric, errors='raise')
        x = numeric.to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError('Numeric raw counts required') from error
    if not np.isfinite(x).all() or (x < 0).any() or not np.equal(x, np.floor(x)).all():
        raise ValueError('Finite nonnegative raw integer counts required')
    # Float conversion must not round larger integers, and group summation must
    # not overflow int64. The bound is conservative for unrealistic huge counts.
    max_count = float(x.max())
    max_cells = int(metadata.groupby('donor', sort=False).size().max())
    if max_count > 2**53 - 1 or int(max_count) * max_cells > np.iinfo(np.int64).max:
        raise ValueError('Counts exceed supported exact integer range')
    if (metadata.groupby('donor', sort=False)[['condition', 'sex']].nunique() > 1).any().any():
        raise ValueError('Donor has conflicting condition or sex')
    numeric = numeric.astype(np.int64)
    grouped = numeric.groupby(metadata['donor'], sort=False).sum()
    donors = metadata.groupby('donor', sort=False)[['condition', 'sex']].first().reindex(grouped.index)
    return grouped, donors


def donor_contrast(counts, metadata, treatment, control, seed=2026, draws=2000):
    """Descriptive natural-log1p count-per-million contrast and donor bootstrap.

    Normalization uses the sum of *supplied* genes. Supply the same full count
    universe for all cells; a filtered panel changes this denominator. This is
    not edgeR, a log2 fold change, an adjusted disease model or a causal estimator.
    Other conditions are rejected rather than silently ignored. Batch, sex,
    disease-stage confounding and repeated/paired donors require another design.
    """
    if treatment == control or isinstance(draws, (bool, np.bool_)) or not isinstance(draws, (int, np.integer)) or draws < 100:
        raise ValueError('Distinct conditions and an integer of at least 100 draws required')
    pb, donors = pseudobulk(counts, metadata)
    if not donors.condition.isin([treatment, control]).all():
        raise ValueError('Subset to the two requested conditions explicitly')
    libraries = pb.astype(float).sum(axis=1)
    if (libraries == 0).any():
        raise ValueError('Zero library-size donor')
    z = np.log1p(pb.div(libraries, axis=0) * 1e6)
    a = z.loc[donors.condition.eq(treatment)].to_numpy()
    b = z.loc[donors.condition.eq(control)].to_numpy()
    if min(len(a), len(b)) < 2:
        raise ValueError('At least two independent donors per group required')
    rng = np.random.default_rng(seed)
    boot = np.array([
        a[rng.integers(len(a), size=len(a))].mean(0)
        - b[rng.integers(len(b), size=len(b))].mean(0)
        for _ in range(draws)
    ])
    return pd.DataFrame({
        'gene': z.columns, 'difference_log1p_cpm': a.mean(0) - b.mean(0),
        'ci_low': np.quantile(boot, .025, axis=0),
        'ci_high': np.quantile(boot, .975, axis=0),
        'n_treatment_donors': len(a), 'n_control_donors': len(b),
    })
