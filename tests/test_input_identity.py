import pandas as pd
import pytest

from athena_atlas.cli import main
from athena_atlas.core import pseudobulk
from athena_atlas.edger import run_edger


def fixture():
    ids = ['001', '01', '002', '02', '003', '03', '004', '04']
    counts = pd.DataFrame({'g': range(1, 9), 'h': [4] * 8}, index=ids)
    meta = pd.DataFrame({'donor': ids, 'condition': ['control'] * 4 + ['treated'] * 4,
                         'sex': ['female'] * 8}, index=ids)
    return counts, meta


def test_csv_keeps_distinct_cell_and_donor_ids(tmp_path):
    counts, meta = fixture()
    counts.to_csv(tmp_path/'counts.csv')
    meta.to_csv(tmp_path/'meta.csv')
    out = tmp_path/'result'
    assert main(['run', '--input', str(tmp_path/'counts.csv'), '--metadata',
                 str(tmp_path/'meta.csv'), '--output', str(out)]) == 0
    donors = pd.read_csv(out/'donors.csv', dtype={'donor': str})
    assert donors.donor.tolist() == meta.donor.tolist()
    contrast = pd.read_csv(out/'contrast.csv')
    assert contrast.n_treatment_donors.tolist() == [4, 4]
    assert contrast.n_control_donors.tolist() == [4, 4]


def test_duplicate_gene_header_cannot_be_renamed(tmp_path):
    counts, meta = fixture()
    counts.columns = ['g', 'g']
    counts.to_csv(tmp_path/'counts.csv')
    meta.to_csv(tmp_path/'meta.csv')
    with pytest.raises(SystemExit) as exc:
        main(['run', '--input', str(tmp_path/'counts.csv'), '--metadata',
              str(tmp_path/'meta.csv'), '--output', str(tmp_path/'result')])
    assert exc.value.code == 2
    assert not (tmp_path/'result').exists()


def test_surrounding_donor_whitespace_rejected():
    counts, meta = fixture()
    meta.iloc[1, 0] = '001 '
    with pytest.raises(ValueError, match='whitespace'):
        pseudobulk(counts, meta)


def test_edger_rejects_repeated_donor_samples_before_r():
    counts, meta = fixture()
    samples = meta.rename(columns={'donor': 'donor_id', 'condition': 'diagnosis',
                                  'sex': 'sex_route'}).assign(sample_id=meta.index)
    samples.iloc[1, 0] = samples.iloc[0, 0]
    with pytest.raises(ValueError, match='one sample per independent donor'):
        run_edger(counts.T, samples, 'female', rscript='unneeded-Rscript')
