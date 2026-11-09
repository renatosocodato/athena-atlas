"""Independent R aggregation and real edgeR interface checks on synthetic counts.

Requires Rscript with edgeR/limma. Missing dependencies are failures of this
explicit integration job, not silently skipped core tests.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from athena_atlas.core import pseudobulk, donor_contrast
from athena_atlas.cli import main as atlas_cli

REFERENCE_R = r'''
args <- commandArgs(trailingOnly=TRUE)
counts <- as.matrix(read.csv(args[1],row.names=1,check.names=FALSE))
meta <- read.csv(args[2],colClasses="character",check.names=FALSE)
rownames(meta) <- meta$cell
stopifnot(identical(rownames(counts),rownames(meta)))
pb <- rowsum(counts,group=meta$donor,reorder=TRUE)
write.csv(pb,args[3])
annotation <- meta[match(rownames(pb),meta$donor),]
norm <- log1p(pb/rowSums(pb)*1e6)
effect <- colMeans(norm[annotation$condition=="treated",,drop=FALSE])-colMeans(norm[annotation$condition=="control",,drop=FALSE])
write.csv(data.frame(gene=names(effect),effect=effect),args[4],row.names=FALSE)
'''

EDGER_REFERENCE_R = r'''
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly=TRUE)
counts <- as.matrix(read.csv(args[1],row.names=1,check.names=FALSE))
m <- read.csv(args[2],stringsAsFactors=FALSE,colClasses=c(sample_id="character",donor_id="character"))
stopifnot(identical(colnames(counts),m$sample_id))
group <- factor(m$diagnosis,levels=c("control","ad"))
y <- DGEList(counts=counts,group=group)
keep <- filterByExpr(y,group=group)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group+braak_numeric,data=m)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
test <- glmQLFTest(fit,coef=2)
tab <- topTags(test,n=Inf,sort.by="none")$table
tab$gene <- rownames(tab)
write.csv(tab,args[3],row.names=FALSE)
writeLines(c(paste("R",getRversion()),paste("edgeR",packageVersion("edgeR")),paste("limma",packageVersion("limma"))),args[4])
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--rscript',default='Rscript')
    args=parser.parse_args()
    if args.output.exists():parser.error('Output directory must be new')
    args.output.mkdir(parents=True)
    output=args.output.resolve()
    rng=np.random.default_rng(1837)
    metadata=[];values=[]
    for donor in range(12):
        for cell in range(3+donor%7):
            metadata.append({'cell':f'cell_{donor}_{cell}','donor':f'{donor:04d}',
                             'condition':'control' if donor<6 else 'treated',
                             'sex':'female' if donor%2==0 else 'male'})
            values.append(rng.negative_binomial(10,.25,size=8))
    meta=pd.DataFrame(metadata).set_index('cell')
    counts=pd.DataFrame(values,index=meta.index,columns=[f'g{i}' for i in range(8)])
    counts.to_csv(output/'cell_counts.csv');meta.reset_index().to_csv(output/'cell_metadata.csv',index=False)
    (output/'aggregation_reference.R').write_text(REFERENCE_R)
    process=subprocess.run([args.rscript,str(output/'aggregation_reference.R'),str(output/'cell_counts.csv'),
                            str(output/'cell_metadata.csv'),str(output/'reference_pb.csv'),str(output/'reference_effect.csv')],capture_output=True,text=True,check=True,timeout=120)
    (output/'aggregation_reference.log').write_text(process.stdout+process.stderr)
    actual,donors=pseudobulk(counts,meta)
    reference=pd.read_csv(output/'reference_pb.csv',dtype={'Unnamed: 0':str}).set_index('Unnamed: 0')
    actual.index=actual.index.astype(str)
    np.testing.assert_array_equal(actual.loc[reference.index].to_numpy(),reference.to_numpy())
    effect=donor_contrast(counts,meta,'treated','control',draws=100)
    expected=pd.read_csv(output/'reference_effect.csv').set_index('gene').effect
    np.testing.assert_allclose(effect.set_index('gene').difference_log1p_cpm.loc[expected.index],expected,rtol=1e-11,atol=1e-11)
    # A separate synthetic RNA-seq fixture exercises real edgeR, not the
    # descriptive Python contrast (the estimators are intentionally different).
    rng=np.random.default_rng(90210)
    means=np.full((200,12),80.)
    means[:20,6:]*=2
    nb=rng.negative_binomial(20,20/(20+means))
    sample_ids=[f'{i:0{width}d}' for i in range(1,7) for width in (3,2)]
    sample_meta=pd.DataFrame({'sample_id':sample_ids,'donor_id':[f'd{i:02d}' for i in range(12)],
                             'sex_route':['female']*12,'diagnosis':['control']*6+['ad']*6,
                             'braak_numeric':[1,2,3,4,5,6]*2})
    frame=pd.DataFrame(nb,index=[f'gene_{i:03d}' for i in range(200)],columns=sample_meta.sample_id)
    frame.to_csv(output/'edger_counts.csv');sample_meta.to_csv(output/'edger_metadata.csv',index=False)
    atlas_cli(['edger','--counts',str(output/'edger_counts.csv'),'--metadata',
               str(output/'edger_metadata.csv'),'--route','female','--rscript',args.rscript,
               '--output',str(output/'adapter')])
    observed=pd.read_csv(output/'adapter/edger.csv')
    (output/'edger_reference.R').write_text(EDGER_REFERENCE_R)
    process=subprocess.run([args.rscript,str(output/'edger_reference.R'),str(output/'edger_counts.csv'),
                            str(output/'edger_metadata.csv'),str(output/'edger_reference.csv'),str(output/'r_versions.txt')],capture_output=True,text=True,check=True,timeout=120)
    (output/'edger_reference.log').write_text(process.stdout+process.stderr)
    expected=pd.read_csv(output/'edger_reference.csv').set_index('gene')
    observed=observed.set_index('gene').loc[expected.index]
    columns=['logFC','logCPM','F','PValue','FDR']
    np.testing.assert_allclose(observed[columns],expected[columns],rtol=1e-10,atol=1e-12)
    summary={'data_kind':'synthetic','donors':len(donors),'cell_counts':len(counts),
             'aggregation':'exact equality with independent R rowsum',
             'descriptive_effect':'matches R natural-log1p CPM mean contrast within 1e-11',
             'edger_genes':len(expected),'edger_adapter_reference':'logFC, logCPM, F, PValue, FDR within rtol 1e-10 / atol 1e-12',
             'r_versions':(output/'r_versions.txt').read_text().splitlines(),
             'interpretation':'Numerical and interface equivalence on synthetic counts; not human atlas validation, sex effect or causal mechanism'}
    (output/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
