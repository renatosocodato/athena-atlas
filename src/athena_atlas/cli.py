DESCRIPTION = 'Donor-aware atlas analysis'

import argparse
import csv
import json
from pathlib import Path
import shutil
import tempfile
import pandas as pd
from . import core
from .edger import run_edger


def _read_csv(path, *, indexed=False):
    """Read identifier fields literally and reject ambiguous headers/rows."""
    with path.open(encoding='utf-8-sig', newline='') as stream:
        rows = csv.reader(stream, strict=True)
        header = next(rows, [])
        names = header[1:] if indexed else header
        if not names or any(not name.strip() for name in names) or len(set(header)) != len(header):
            raise ValueError('CSV requires nonempty unique column names')
        records = []
        for row in rows:
            if not row:
                continue
            if len(row) != len(header) or (indexed and not row[0].strip()):
                raise ValueError(f'Invalid CSV record at line {rows.line_num}')
            records.append(row)
    frame = pd.DataFrame(records, columns=header)
    return frame.set_index(header[0]) if indexed else frame

def main(argv=None):
    parser=argparse.ArgumentParser(description=DESCRIPTION)
    sub=parser.add_subparsers(dest='command',required=True)
    ex=sub.add_parser('example',help='Run a synthetic example, without network or private data')
    ex.add_argument('--output',type=Path,required=True)
    run=sub.add_parser('run',help='Analyze explicit input files; see INPUTS.md')
    run.add_argument('--input',type=Path,required=True);run.add_argument('--output',type=Path,required=True)
    run.add_argument('--metadata',type=Path);run.add_argument('--treatment',default='treated');run.add_argument('--control',default='control')
    edger = sub.add_parser('edger', help='Run edgeR on one raw pseudobulk sample per independent donor')
    edger.add_argument('--counts', type=Path, required=True)
    edger.add_argument('--metadata', type=Path, required=True)
    edger.add_argument('--route', required=True)
    edger.add_argument('--output', type=Path, required=True)
    edger.add_argument('--rscript', default='Rscript')
    args=parser.parse_args(argv)
    if args.output.exists():parser.error('Output directory already exists; choose a new directory')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(prefix='.athena-stage-',dir=args.output.parent) as tmp:
            staged=argparse.Namespace(**vars(args));staged.output=Path(tmp)
            result=analyze(staged)
            (staged.output/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
            shutil.copytree(staged.output,args.output)
    except (ValueError,KeyError,OSError,RuntimeError,csv.Error) as error:
        parser.error(str(error))
    return 0

def analyze(args):
    if args.command == 'edger':
        frame = run_edger(_read_csv(args.counts, indexed=True), _read_csv(args.metadata),
                          args.route, rscript=args.rscript)
        frame.to_csv(args.output/'edger.csv', index=False)
        return {'data_kind': 'user_supplied', 'genes': len(frame),
                'method': 'edgeR quasi-likelihood donor pseudobulk contrast'}
    if args.command=='example':
        metadata=pd.DataFrame({'donor':[f'd{i}' for i in range(8) for _ in range(3)],'condition':['control']*12+['treated']*12,'sex':['female']*24},index=[f'cell{i}' for i in range(24)])
        counts=pd.DataFrame({'gene_a':[5]*12+[12]*12,'gene_b':[12]*12+[5]*12},index=metadata.index)
    else:
        if args.metadata is None:raise ValueError('--metadata required')
        counts=_read_csv(args.input,indexed=True);metadata=_read_csv(args.metadata,indexed=True)
    pb,donors=core.pseudobulk(counts,metadata)
    contrast=core.donor_contrast(counts,metadata,args.treatment if args.command=='run' else 'treated',args.control if args.command=='run' else 'control')
    pb.to_csv(args.output/'pseudobulk.csv');donors.to_csv(args.output/'donors.csv');contrast.to_csv(args.output/'contrast.csv',index=False)
    return {'data_kind':'synthetic' if args.command=='example' else 'user_supplied','donors':len(donors),'method':'descriptive log1p-CPM contrast; donor bootstrap; no covariate adjustment'}
