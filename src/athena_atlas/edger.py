"""Validated adapter to an explicitly installed R/edgeR backend."""
from importlib.resources import as_file, files
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd


def run_edger(counts, metadata, route, *, rscript='Rscript'):
    """Fit an AD/control contrast to gene-by-donor raw counts in one route.

    Each sample must represent one independent donor. The selected route must
    contain both diagnoses and may include a finite Braak covariate. This adapter
    validates correspondence and leaves biological design choices to the caller.
    """
    if not isinstance(counts, pd.DataFrame) or not isinstance(metadata, pd.DataFrame):
        raise ValueError('Counts and metadata must be data frames')
    required = ['sample_id', 'donor_id', 'sex_route', 'diagnosis']
    if (counts.empty or not counts.index.is_unique or not counts.columns.is_unique
            or counts.index.hasnans or counts.columns.hasnans or not metadata.columns.is_unique
            or not set(required).issubset(metadata.columns)):
        raise ValueError('Unique gene/sample IDs and complete donor metadata required')
    metadata = metadata.copy()
    for column in required:
        values = metadata[column].astype('string')
        if values.isna().any() or values.str.strip().eq('').any():
            raise ValueError(f'Nonempty {column} required')
        metadata[column] = values
    genes = counts.index.map(str)
    samples = counts.columns.map(str)
    if (not genes.is_unique or not samples.is_unique or any(not value.strip() for value in genes)
            or metadata.sample_id.duplicated().any() or metadata.donor_id.duplicated().any()
            or set(samples) != set(metadata.sample_id)):
        raise ValueError('Counts must match one sample per independent donor')
    if not isinstance(route, str) or not route.strip():
        raise ValueError('A nonempty route is required')
    selected = metadata.set_index('sample_id').loc[samples].reset_index(names='sample_id')
    selected = selected.loc[selected.sex_route.eq(route)].copy()
    diagnoses = selected.diagnosis.value_counts().to_dict()
    if set(diagnoses) != {'control', 'ad'} or min(diagnoses.values()) < 2:
        raise ValueError('Selected route requires at least two independent donors per diagnosis')
    if 'braak_numeric' in selected:
        selected['braak_numeric'] = pd.to_numeric(selected.braak_numeric, errors='raise')
        if not np.isfinite(selected.braak_numeric.to_numpy(dtype=float)).all():
            raise ValueError('Braak covariate must be finite when supplied')
    frame = counts.copy(deep=False)
    frame.index, frame.columns = genes, samples
    numeric = frame.loc[:, selected.sample_id].apply(pd.to_numeric, errors='raise')
    values = numeric.to_numpy(dtype=float)
    if (not np.isfinite(values).all() or (values < 0).any() or (values > 2**53-1).any()
            or not np.equal(values, np.floor(values)).all() or (values.sum(axis=0) == 0).any()):
        raise ValueError('Finite nonnegative integer raw counts and positive libraries required')
    with tempfile.TemporaryDirectory(prefix='athena_edger_') as directory:
        folder = Path(directory)
        counts_path, metadata_path, result = (folder / name for name in ('counts.csv', 'metadata.csv', 'result.csv'))
        numeric.astype(np.int64).to_csv(counts_path)
        selected.to_csv(metadata_path, index=False)
        with as_file(files('athena_atlas').joinpath('resources', 'edger_route.R')) as script:
            try:
                subprocess.run([rscript, str(script), str(counts_path), str(metadata_path), route, str(result)],
                               check=True, capture_output=True, text=True, timeout=120)
            except subprocess.CalledProcessError as error:
                raise RuntimeError(error.stderr.strip() or 'edgeR failed') from error
            except subprocess.TimeoutExpired as error:
                raise RuntimeError('edgeR exceeded its 120 second timeout') from error
        if not result.is_file():
            raise RuntimeError('edgeR did not produce a result')
        output = pd.read_csv(result, dtype={'gene': str, 'sex_route': str}, keep_default_na=False)
        if output.empty:
            raise RuntimeError('edgeR produced an empty result')
        return output
