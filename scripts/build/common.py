"""Shared helpers for the database builders.

Every builder produces a fresh SQLite file, stamps it with a one-row
`_version` provenance table, and writes tables whose schemas and naming
conventions match what curie's readers require (see the curie repository's
`Library`/`Isotope` classes for the consuming SQL).
"""
import os
import re
import sqlite3
import subprocess
from datetime import date
from pathlib import Path

import numpy as np

ARCHIVE_ROOT = Path(os.environ.get(
    'CURIE_SOURCE_ARCHIVE',
    Path.home() / 'Documents/nucleardata/data_sources/v2'))

REPO_ROOT = Path(__file__).resolve().parents[2]

DECAY_DB = Path(os.environ.get('CURIE_DECAY_DB', REPO_ROOT / 'build' / 'decay.db'))


def curated_targets():
    """The curated TENDL target set (curie-style names) from tendl_targets.txt."""
    path = Path(__file__).parent / 'tendl_targets.txt'
    return {t for t in path.read_text().split() if t}


def isomer_ordinal(letter):
    """TALYS isomer-target letter -> curie ordinal: '' -> 0, 'm' -> 1,
    'n' -> 2, 'o' -> 3 (the inverse lives in scripts/fetch/tendl.py's
    tendl_paths, which builds TENDL paths from curie names)."""
    return ord(letter) - ord('m') + 1 if letter else 0


def ordinal_suffix(k):
    """curie isomer-target suffix for ordinal k: 0 -> '', 1 -> 'm' (the
    established bare-m first-isomer convention), 2 -> 'm2', ..."""
    return '' if not k else 'm' if k == 1 else 'm{}'.format(k)


def parse_mirror_dir(dirname):
    """'Hf178n' -> ('Hf', 178, 'm2'): element, mass number, and curie
    isomer suffix of a TENDL mirror directory name."""
    m = re.fullmatch(r'([A-Z][a-z]?)(\d{3})([mno]?)', dirname)
    return m.group(1), int(m.group(2)), ordinal_suffix(isomer_ordinal(m.group(3)))


def mirror_target(dirname):
    """A mirror directory name as a curie-style target: 'Ag108m' ->
    '108AGm', 'Hf178n' -> '178HFm2'."""
    el, a, iso = parse_mirror_dir(dirname)
    return f'{a}{el.upper()}{iso}'


BUILD_DIR = REPO_ROOT / 'build'


def open_fresh(name):
    """Create build/<name> anew and return a connection."""
    BUILD_DIR.mkdir(exist_ok=True)
    path = BUILD_DIR / name
    if path.exists():
        path.unlink()
    return sqlite3.connect(path)


def pipeline_commit():
    try:
        out = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
                             cwd=REPO_ROOT, capture_output=True, text=True)
        sha = out.stdout.strip()
        dirty = subprocess.run(['git', 'status', '--porcelain'],
                               cwd=REPO_ROOT, capture_output=True, text=True)
        return sha + ('+dirty' if dirty.stdout.strip() else '')
    except Exception:
        return 'unknown'


def write_version(con, library, source, source_version, generation='v2'):
    """One-row provenance table recording what this file was built from.

    Nothing reads it programmatically yet; it exists so the file itself
    carries its generation and source (inspectable with any sqlite tool),
    and so a future curie release can detect generation skew.
    """
    con.execute('CREATE TABLE _version (library TEXT, generation TEXT, '
                'source TEXT, source_version TEXT, build_date TEXT, '
                'pipeline TEXT)')
    con.execute('INSERT INTO _version VALUES (?,?,?,?,?,?)',
                (library, generation, source, source_version,
                 date.today().isoformat(), pipeline_commit()))


# --- nuclide / reaction-string normalization (curie conventions) ---

# Same ordering as curie's Isotope class: list index = Z.
ELEMENTS = ['n','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
            'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',
            'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr',
            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I',
            'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
            'Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au',
            'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U',
            'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db',
            'Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts']

def curie_nuclide(token):
    """Nuclide token -> curie name, accepting both source spellings:

    mass-first (IAEA):    '22Na' -> '22NA', '94mTc' -> '94TCm',
                          '102gRh' -> '102RHg', '52Mnm' -> '52MNm'
    element-first (TALYS residual headers): 'Mn52m' -> '52MNm',
                          'Co57' -> '57CO'

    Tokens with no isomer marker keep no suffix; two markers at once are
    rejected.
    """
    import re
    t = token.strip()
    m = re.fullmatch(r'(\d+)(m\d*|g)?([A-Za-z]{1,2})(m\d*|g)?', t)
    if m:
        a, iso1, el, iso2 = m.groups()
    else:
        m = re.fullmatch(r'([A-Za-z]{1,2})(\d+)(m\d*|g)?', t)
        if not m:
            raise ValueError(f'unparseable nuclide token: {token!r}')
        el, a, iso1 = m.groups()
        iso2 = None
    if iso1 and iso2:
        raise ValueError(f'two isomer markers in nuclide token: {token!r}')
    return f'{a}{el.upper()}{iso1 or iso2 or ""}'


def parse_yandf(text):
    """(header_dict, rows[(E [MeV], sigma [mb])]) from a TALYS YANDF table.

    Header lines are '#   key: value' pairs (nested sections flattened,
    last occurrence wins); data rows take the first two columns.
    """
    header, rows = {}, []
    for line in text.splitlines():
        if line.startswith('#'):
            if ':' in line:
                k, _, v = line.lstrip('# ').partition(':')
                header[k.strip()] = v.strip()
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    return header, rows


def dedupe_grid(rows, where, rel_tol=1e-5):
    """Collapse duplicate abscissae whose values agree to `rel_tol`
    (round-off-level duplicated print rows are endemic in the TALYS
    tables); a duplicate disagreeing beyond tolerance is a data conflict
    and fails the build. Returns sorted (x, y) pairs."""
    seen = {}
    for e, x in rows:
        if e in seen:
            ref = seen[e]
            if abs(x - ref) > rel_tol * max(abs(ref), abs(x), 1e-300):
                raise SystemExit(f'conflicting duplicate abscissa {e} in '
                                 f'{where}: {ref} vs {x}')
            continue  # keep-first among round-off twins
        seen[e] = x
    return sorted(seen.items())


def union_grid_stack(series):
    """Assemble channel/product series onto their union energy grid.

    `series` is a list of rows[(E, value)]; returns (grid, data) where data
    is the column stack [grid, col0, col1, ...] with each column linearly
    interpolated onto the union grid and zero-filled outside its own
    tabulated range (the TALYS production files are threshold-truncated
    subsets of the master grid).
    """
    grid = np.unique(np.concatenate([[e for e, _ in r] for r in series]))
    data = np.column_stack(
        [grid] + [np.interp(grid, [e for e, _ in r], [x for _, x in r],
                            left=0.0, right=0.0) for r in series])
    return grid, data


def curie_target(token):
    """'natTi' -> 'natTI', '27Al' -> '27AL', '60Ni' -> '60NI'."""
    token = token.strip()
    if token.lower().startswith('nat'):
        return 'nat' + token[3:].upper()
    return curie_nuclide(token)
