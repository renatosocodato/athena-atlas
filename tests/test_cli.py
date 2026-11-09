import subprocess
import sys
from pathlib import Path
import pytest

def test_example_and_no_overwrite(tmp_path):
    out=tmp_path/'example'
    command=[sys.executable,'-m','athena_atlas','example','--output',str(out)]
    subprocess.run(command,check=True,capture_output=True)
    assert (out/'result.json').is_file()
    second=subprocess.run(command,capture_output=True)
    assert second.returncode!=0


def test_invalid_input_leaves_no_output_or_staging_directory(tmp_path):
    bad=tmp_path/'bad.csv';bad.write_text('wrong,column\n1,2\n')
    out=tmp_path/'out'
    command=[sys.executable,'-m','athena_atlas','run','--input',str(bad),'--output',str(out)]
    result=subprocess.run(command,capture_output=True)
    assert result.returncode!=0
    assert not out.exists()
    assert not list(tmp_path.glob('.athena-stage-*'))
