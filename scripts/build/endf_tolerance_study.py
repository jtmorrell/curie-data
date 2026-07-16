"""RECONR linearization-tolerance study for the endf.db rebuild.

Runs NJOY RECONR on a set of resonance-heavy materials at several
fractional tolerances, reports pointwise grid sizes, and projects the
whole-library database size from the previous generation's row counts.
The result picks the production tolerance.

Usage: python3 endf_tolerance_study.py
Requires: the ENDF/B-VIII.1 neutrons tarball in the source archive and
`njoy` on PATH or at ~/micromamba/envs/njoy/bin/njoy.
"""
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import ARCHIVE_ROOT

NJOY = shutil.which('njoy') or str(Path.home() / 'micromamba/envs/njoy/bin/njoy')
TARBALL = ARCHIVE_ROOT / 'endf-b-viii.1' / 'neutrons-version.VIII.1.tar.gz'
V1_DIR = Path(os.environ.get('CURIE_V1_DIR', Path.home() / '.local/share/curie'))

TARGETS = ['n-026_Fe_056.endf', 'n-048_Cd_113.endf', 'n-092_U_235.endf',
           'n-074_W_186.endf', 'n-079_Au_197.endf']
TOLERANCES = [0.001, 0.005, 0.01]


def mat_number(endf_path):
    """MAT from the first MF1/MT451 card (columns 67-70)."""
    with open(endf_path) as f:
        f.readline()  # tape ID
        return int(f.readline()[66:70])


def reconr_points(endf_path, tol, workdir):
    """Total pointwise rows RECONR produces at `tol` (sum over MF3 MTs)."""
    shutil.copy(endf_path, workdir / 'tape20')
    mat = mat_number(endf_path)
    deck = (f"reconr\n20 21\n'tolerance study'/\n{mat} 1/\n{tol}/\n"
            f"'{endf_path.name}'/\n0/\nstop\n")
    subprocess.run([NJOY], input=deck, text=True, cwd=workdir,
                   capture_output=True, check=True)
    import endf
    mat_out = endf.Material(str(workdir / 'tape21'))
    energies = mat_out.section_data[(3, 1)]['sigma'].x
    total = sum(len(s['sigma'].x) for (mf, mt), s in mat_out.section_data.items()
                if mf == 3 and 'sigma' in s)
    (workdir / 'tape21').unlink()
    return len(energies), total


def main():
    with tarfile.open(TARBALL) as t:
        members = {Path(m.name).name: m for m in t.getmembers()}
        tmp = Path(tempfile.mkdtemp(prefix='reconr-study-'))
        for name in TARGETS:
            t.extract(members[name], tmp, filter='data')

    v1 = sqlite3.connect(V1_DIR / 'endf.db')
    v1_rows = {t: v1.execute(f'SELECT COUNT(*) FROM {t}').fetchone()[0]
               for t in ['FE_56', 'CD_113', 'U_235', 'W_186', 'AU_197']}
    v1_total = sum(r[0] for r in v1.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name!='all_reactions'")
        for r in [v1.execute(f'SELECT COUNT(*) FROM {r[0]}').fetchone()])

    print(f'{"material":<12}{"v1 rows":>9}' +
          ''.join(f'{"tol=" + str(t):>14}' for t in TOLERANCES))
    ratios = {t: [] for t in TOLERANCES}
    for name in TARGETS:
        path = next(tmp.rglob(name))
        key = name.split('_', 1)[1].replace('.endf', '').replace('_0', '_') \
                  .replace('Fe_', 'FE_').upper()
        line = f'{key:<12}{v1_rows.get(key, 0):>9}'
        for tol in TOLERANCES:
            with tempfile.TemporaryDirectory() as wd:
                n_e, _ = reconr_points(path, tol, Path(wd))
            line += f'{n_e:>14}'
            if v1_rows.get(key):
                ratios[tol].append(n_e / v1_rows[key])
        print(line, flush=True)

    print(f'\nv1 endf.db total rows: {v1_total} (748 MB whole file)')
    for tol in TOLERANCES:
        r = sum(ratios[tol]) / len(ratios[tol])
        print(f'tol={tol}: mean grid ratio vs v1 = {r:.2f} '
              f'-> projected whole file ~{int(748 * r)} MB')


if __name__ == '__main__':
    main()
