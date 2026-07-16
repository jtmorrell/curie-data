"""Shared helpers for the database builders.

Every builder produces a fresh SQLite file, stamps it with a one-row
`_version` provenance table, and writes tables whose schemas and naming
conventions match what curie's readers require (see the curie repository's
`Library`/`Isotope` classes for the consuming SQL).
"""
import os
import sqlite3
import subprocess
from datetime import date
from pathlib import Path

ARCHIVE_ROOT = Path(os.environ.get(
    'CURIE_SOURCE_ARCHIVE',
    Path.home() / 'Documents/nucleardata/data_sources/v2'))

REPO_ROOT = Path(__file__).resolve().parents[2]
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


def curie_target(token):
    """'natTi' -> 'natTI', '27Al' -> '27AL', '60Ni' -> '60NI'."""
    token = token.strip()
    if token.lower().startswith('nat'):
        return 'nat' + token[3:].upper()
    return curie_nuclide(token)
