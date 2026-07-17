"""Build endf.db from the ENDF/B-VIII.1 neutron sublibrary.

Per material: RECONR reconstructs the resonance ranges to 0 K pointwise
(fractional tolerance below), an optional HEATR pass adds the MT444 damage
section, and the MF3 sections are assembled into one wide table on the
union energy grid (the MT1 grid): one `energy` column [eV] plus one column
per stored channel [barn], zero-filled outside a channel's range. Channel
selection, column labels, outgoing strings, and residual products come
from endf_mt_labels.toml; the union of MF3 sections neither mapped nor
skip-listed is reported at the end of the run.

The 'abs' column is computed as sum(MT102..117) — the MT101 disappearance
definition, excluding fission, matching the previous generation. A HEATR
failure drops only that material's 'damage' column (counted in the report).

The build is resumable: each material commits atomically (table + rows +
index entries), and a material whose table exists is skipped, so an
interrupted run continues where it left off; the end-of-run accounting
covers the current invocation only. Resuming with a different --tolerance
is refused. Table names derive from the source file names
('n-048_Cd_113.endf' -> CD_113, 'n-095_Am_242m1.endf' -> AM_242m).

Usage: python3 endf_neutrons.py [--tolerance 0.005] [--limit N]
"""
import re
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
import tomllib
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from common import ARCHIVE_ROOT, BUILD_DIR, ELEMENTS, write_version

NJOY = shutil.which('njoy') or str(Path.home() / 'micromamba/envs/njoy/bin/njoy')
TARBALL = ARCHIVE_ROOT / 'endf-b-viii.1' / 'neutrons-version.VIII.1.tar.gz'
MAPPING = tomllib.load(open(Path(__file__).parent / 'endf_mt_labels.toml', 'rb'))

PARTICLE_AZ = {'n': (1, 0), 'p': (1, 1), 'd': (2, 1), 't': (3, 1),
               'h': (3, 2), 'a': (4, 2)}
# absorption = MT101 disappearance (sum of MT102..117), matching the
# previous generation; fission is NOT included (it lives in _fis_)
ABS_MTS = list(range(102, 118))


def mt_plan():
    """{mt: (label, outgoing, kind)}; kind 'summary'|'emission'|'level'|'production'."""
    plan = {}
    for lab, mt in MAPPING['summary'].items():
        out = {'fis': 'f', 'f': 'f1'}.get(lab, lab)
        plan[mt] = (f'_{lab}_', out, 'summary')
    for lab, spec in MAPPING['emission'].items():
        plan[spec['mt']] = (f'_{lab}_', lab, 'emission')
    for fam, spec in MAPPING['levels'].items():
        # the n-family level 0 (MT50) does not occur; the charged-particle
        # families' level 0 (MT600/650/700/750/800 -> _p0_ etc.) is real
        start = 1 if fam == 'n' else 0
        for i in range(start, spec['count']):
            plan[spec['first'] + i] = (f'_{fam}{i}_', f'{fam}{i}', 'level')
        plan[spec['continuum']] = (f'_{fam}c_', f'{fam}c', 'level')
    for lab, mt in MAPPING['production'].items():
        plan[mt] = (f'_{lab}', 'x', 'production')
    return plan


MT_PLAN = mt_plan()
SKIP = set(MAPPING['skip']['mts'])


def material_name(fname):
    """'n-048_Cd_113.endf' -> ('CD_113', A=113, Z=48); m1 -> trailing 'm'."""
    m = re.fullmatch(r'n-(\d+)_([A-Za-z]+)_(\d+)(m\d*)?\.endf', fname)
    z, el, a, iso = int(m.group(1)), m.group(2), int(m.group(3)), m.group(4)
    return f'{el.upper()}_{a}{"m" if iso else ""}', a, z


def residual(a, z, emit):
    """Ground-state product name for target (a, z) + n - emitted particles."""
    pa, pz = a + 1, z
    for p in emit.split():
        da, dz = PARTICLE_AZ[p]
        pa, pz = pa - da, pz - dz
    return f'{pa}{ELEMENTS[pz].upper()}' if 0 < pz < len(ELEMENTS) else ''


def njoy_pointwise(endf_path, tol, workdir, want_damage=True):
    """RECONR (+HEATR) -> parsed PENDF Material; damage flag in second slot."""
    shutil.copy(endf_path, workdir / 'tape20')
    with open(endf_path) as f:
        f.readline()
        mat = int(f.readline()[66:70])
    deck = (f"reconr\n20 21\n'pointwise'/\n{mat} 1/\n{tol}/\n''/\n0/\nstop\n")
    subprocess.run([NJOY], input=deck, text=True, cwd=workdir,
                   capture_output=True, check=True)
    tape = 'tape21'
    damage = False
    if want_damage:
        deck = f"heatr\n20 21 22/\n{mat} 1/\n444/\nstop\n"
        r = subprocess.run([NJOY], input=deck, text=True, cwd=workdir,
                           capture_output=True)
        if r.returncode == 0 and (workdir / 'tape22').exists():
            tape, damage = 'tape22', True
    import endf
    return endf.Material(str(workdir / tape)), damage


def build_material(con, endf_path, tol):
    table, a, z = material_name(endf_path.name)
    with tempfile.TemporaryDirectory() as wd:
        mat, has_damage = njoy_pointwise(endf_path, tol, Path(wd))
    sections = {mt: s['sigma'] for (mf, mt), s in mat.section_data.items()
                if mf == 3 and 'sigma' in s}
    # RECONR occasionally emits duplicated grid points (identical rows, not
    # step discontinuities); np.unique sorts and collapses them
    grid = np.unique(sections[1].x)
    cols, unmapped = [], []
    for mt in sorted(sections):
        if mt in SKIP or mt == 444 and not has_damage:
            continue
        if mt not in MT_PLAN:
            unmapped.append(mt)
            continue
        label, out, kind = MT_PLAN[mt]
        sig = sections[mt]
        y = np.interp(grid, sig.x, sig.y, left=0.0, right=0.0)
        if kind == 'emission':
            prod = residual(a, z, MAPPING['emission'][out]['emit'])
        elif kind == 'production':
            prod = label.rsplit('_', 1)[-1]
        else:
            prod = ''
        cols.append((label, out, prod, y))
    # computed absorption on the union grid
    present = [m for m in ABS_MTS if m in sections]
    if present:
        y = sum(np.interp(grid, sections[m].x, sections[m].y, left=0.0,
                          right=0.0) for m in present)
        cols.insert(1, ('_abs_', 'abs', '', y))

    # one transaction per material: the table, its rows, and its index
    # entries land together or not at all. The explicit BEGIN IMMEDIATE is
    # load-bearing — CREATE TABLE is DDL, which sqlite3 autocommits outside
    # an open transaction, so without it a killed run leaves a hollow table
    # for the resume check to mistake for done work.
    col_sql = ', '.join(f'"{lab}" REAL' for lab, *_ in cols)
    try:
        con.execute('BEGIN IMMEDIATE')
        con.execute(f'CREATE TABLE {table} (energy REAL, {col_sql})')
        data = np.column_stack([grid] + [y for *_, y in cols])
        con.executemany(
            f'INSERT INTO {table} VALUES ({",".join("?" * (len(cols) + 1))})',
            data.tolist())
        tgt = f'{a}{ELEMENTS[z].upper()}' + ('m' if table.endswith('m') else '')
        con.executemany('INSERT INTO all_reactions VALUES (?,?,?,?)',
                        [(tgt, out, prod, lab) for lab, out, prod, _ in cols])
        con.commit()
    except BaseException:
        con.rollback()
        raise
    return table, len(grid), len(cols), unmapped, has_damage


def parse_args(argv):
    tol, limit = 0.005, None
    it = iter(argv)
    for a in it:
        if a == '--tolerance':
            tol = float(next(it, 'missing'))
        elif a == '--limit':
            limit = int(next(it, 'missing'))
        else:
            sys.exit(f'unknown argument {a!r}\n\n{__doc__.strip()}')
    return tol, limit


def main():
    try:
        tol, limit = parse_args(sys.argv[1:])
    except ValueError as e:
        sys.exit(f'bad argument value: {e}')

    BUILD_DIR.mkdir(exist_ok=True)
    con = sqlite3.connect(BUILD_DIR / 'endf.db')
    # a rowless material table is a hollow remnant of an interrupted run
    # (possible in files written before the BEGIN IMMEDIATE guard):
    # rebuild it rather than skipping it as done
    have = set()
    for (name,) in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'").fetchall():
        if (name in ('all_reactions', '_build_meta', '_version')
                or con.execute(f'SELECT 1 FROM {name} LIMIT 1').fetchone()):
            have.add(name)
        else:
            con.execute(f'DROP TABLE {name}')
    con.commit()
    if 'all_reactions' not in have:
        con.execute('CREATE TABLE all_reactions (target TEXT, outgoing TEXT, '
                    'product TEXT, label TEXT)')
        con.execute('CREATE TABLE _build_meta (tolerance REAL)')
        con.execute('INSERT INTO _build_meta VALUES (?)', (tol,))
        con.commit()
    else:
        rec = con.execute('SELECT tolerance FROM _build_meta').fetchone()[0]
        if rec != tol:
            sys.exit(f'existing build used tolerance {rec}; resuming at {tol} '
                     f'would mix grids — delete build/endf.db or match it')
    all_unmapped, no_damage, failures = {}, [], []
    n_done = 0
    with tarfile.open(TARBALL) as t:
        members = [m for m in t.getmembers() if m.name.endswith('.endf')]
        print(f'{len(members)} materials; tolerance {tol}')
        for m in sorted(members, key=lambda m: m.name):
            fname = Path(m.name).name
            table = material_name(fname)[0]
            if table in have:
                continue
            with tempfile.TemporaryDirectory() as xd:
                t.extract(m, xd, filter='data')
                src = next(Path(xd).rglob(fname))
                try:
                    table, npts, ncols, unmapped, dmg = build_material(con, src, tol)
                except Exception as e:
                    failures.append((fname, repr(e)))
                    print(f'{fname}: FAILED {e}', flush=True)
                    continue
            if unmapped:
                all_unmapped[table] = unmapped
            if not dmg:
                no_damage.append(table)
            n_done += 1
            if n_done % 25 == 0:
                print(f'... {n_done} built (latest {table}: {npts} pts, '
                      f'{ncols} channels)', flush=True)
            if limit and n_done >= limit:
                break

    tables = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    complete = not failures and not limit and \
        len(tables - {'all_reactions', '_version', '_build_meta'}) == len(members)
    if complete and '_version' not in tables:
        write_version(con, 'endf',
                      'ENDF/B-VIII.1 neutron sublibrary (NNDC), RECONR '
                      f'tolerance {tol}', 'VIII.1 (2024-08-30)')
        con.commit()
    print(f'\ndone: {n_done} built, {len(failures)} failures, '
          f'{len(no_damage)} without damage (HEATR)')
    for f in failures:
        print('  FAILED:', *f)
    if all_unmapped:
        mts = sorted({mt for v in all_unmapped.values() for mt in v})
        print(f'unmapped MTs across {len(all_unmapped)} materials: {mts}')
    con.close()


if __name__ == '__main__':
    main()
