from pathlib import Path
import os
import re
import shlex
import subprocess
import ncpus

from jug.utils import cached_glob
from jug import TaskGenerator


INPUT_DIR = Path("input")
OUTPUT_DIR = Path("output") / "eggnog"
EGGNOG_MAPPER_DB = Path('~/databases/emapper_db').expanduser()




@TaskGenerator
def run_eggnog_mapper(input_file):
    import tempfile
    import shutil
    input_path = Path(input_file)
    prefix = input_path.with_suffix('').stem
    output_dir = OUTPUT_DIR / prefix
    output_dir.mkdir(parents=True, exist_ok=True)

    extra_args = shlex.split(os.environ.get("EGGNOG_MAPPER_ARGS", ""))
    with tempfile.TemporaryDirectory() as tdir:
        edir = tdir + '/emapper_db'
        shutil.copytree(EGGNOG_MAPPER_DB, edir)

        cmd = [
            "emapper.py",
            "-i", str(input_path),
            "--itype", "genome",
            "--genepred", "prodigal",
            "--output", prefix,
            "--output_dir", str(output_dir),
            "--cpu", str(ncpus.get_ncpus()),
            "--data_dir", edir,
            "--temp_dir", tdir,
            "--override",
        ] + extra_args
        print(cmd)

        subprocess.run(cmd, check=True)

    annotations = output_dir / f"{prefix}.emapper.annotations"
    if not annotations.exists():
        raise FileNotFoundError(f"eggNOG-mapper did not create {annotations}")
    return str(annotations)

@TaskGenerator
def create_magsview_data(annot_file):
    annot_file = annot_file.replace('.fna', '')
    import magsviewdata

    summary = magsviewdata.create_emapper_summary(Path(annot_file).parent)
    subprocess.check_call(
            ['xz',
             '--threads', '2',
             annot_file])
    return summary


for p in cached_glob(f'input/*.fna.gz'):
    emap = run_eggnog_mapper(p)
    if 'MAGSVIEW' in os.environ:
        create_magsview_data(emap)
