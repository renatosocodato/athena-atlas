import numpy as np
import pandas as pd
import pytest
from athena_atlas.core import pseudobulk, donor_contrast

def fixture():
    meta=pd.DataFrame({'donor':['a','a','b','b','c','c','d','d'],'condition':['control']*4+['treated']*4,'sex':['female']*8})
    counts=pd.DataFrame({'x':[1,2,3,4,7,8,9,10],'y':[10,9,8,7,4,3,2,1]})
    return counts,meta

def test_exact_pseudobulk_and_independent_donor_n():
    x,m=fixture();pb,donors=pseudobulk(x,m)
    assert pb.loc['a','x']==3
    assert pb.loc['c','y']==7
    result=donor_contrast(x,m,'treated','control',draws=100)
    assert result.n_treatment_donors.tolist()==[2,2]
    assert result.iloc[0].difference_log1p_cpm>0

def test_cell_duplication_does_not_inflate_donor_n_or_change_normalized_contrast():
    x,m=fixture();a=donor_contrast(x,m,'treated','control',draws=100)
    xx=pd.concat([x,x],ignore_index=True);mm=pd.concat([m,m],ignore_index=True)
    b=donor_contrast(xx,mm,'treated','control',draws=100)
    np.testing.assert_allclose(a.difference_log1p_cpm,b.difference_log1p_cpm)
    assert b.n_treatment_donors.tolist()==[2,2]

@pytest.mark.parametrize('invalid',['fraction','negative','conflict','unaligned','zero'])
def test_invalid_inputs(invalid):
    x,m=fixture()
    if invalid=='fraction':x=x.astype(float);x.loc[0,'x']=.5
    if invalid=='negative':x.loc[0,'x']=-1
    if invalid=='conflict':m.loc[0,'sex']='male'
    if invalid=='unaligned':m.index=m.index+10
    if invalid=='zero':x.loc[:1,:]=0
    with pytest.raises(ValueError):donor_contrast(x,m,'treated','control',draws=100)


def test_numeric_strings_sum_as_counts_not_concatenated_text():
    x,m=fixture();pb,_=pseudobulk(x.astype(str),m)
    assert pb.loc['a','x']==3
    assert pb.loc['d','x']==19

@pytest.mark.parametrize('invalid',['empty','duplicate_genes','missing_gene','empty_gene','missing_metadata','blank_donor','missing_cell','huge_counts'])
def test_schema_failures_are_explicit(invalid):
    x,m=fixture()
    if invalid=='empty':x=x.iloc[:0];m=m.iloc[:0]
    if invalid=='duplicate_genes':x.columns=['x','x']
    if invalid=='missing_gene':x.columns=['x',None]
    if invalid=='empty_gene':x.columns=['x','  ']
    if invalid=='missing_metadata':m=m.drop(columns=['sex'])
    if invalid=='blank_donor':m.loc[0,'donor']='  '
    if invalid=='missing_cell':x.index=[None,*range(1,len(x))];m.index=x.index
    if invalid=='huge_counts':x.loc[0,'x']=2**53
    with pytest.raises(ValueError):pseudobulk(x,m)

@pytest.mark.parametrize('draws',[100.5,True,99])
def test_bootstrap_requires_valid_integer_draw_count(draws):
    x,m=fixture()
    with pytest.raises(ValueError):donor_contrast(x,m,'treated','control',draws=draws)

def test_extra_condition_not_silently_discarded():
    x,m=fixture();m.loc[m.donor.eq('b'),'condition']='other'
    with pytest.raises(ValueError,match='two requested'):donor_contrast(x,m,'treated','control',draws=100)
