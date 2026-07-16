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
    """One-row provenance table; curie reads it at connection open."""
    con.execute('CREATE TABLE _version (library TEXT, generation TEXT, '
                'source TEXT, source_version TEXT, build_date TEXT, '
                'pipeline TEXT)')
    con.execute('INSERT INTO _version VALUES (?,?,?,?,?,?)',
                (library, generation, source, source_version,
                 date.today().isoformat(), pipeline_commit()))


# --- nuclide / reaction-string normalization (curie conventions) ---

def curie_nuclide(token):
    """IAEA/TALYS nuclide token -> curie name: '22Na' -> '22NA',
    '94mTc' -> '94TCm', '102gRh' -> '102RHg', 'Mn52m' -> '52MNm'."""
    import re
    m = re.fullmatch(r'(\d+)(m\d*|g)?([A-Za-z]{1,3})(m\d*|g)?', token.strip())
    if not m:
        raise ValueError(f'unparseable nuclide token: {token!r}')
    a, iso1, el, iso2 = m.groups()
    iso = iso1 or iso2 or ''
    return f'{a}{el.upper()}{iso}'


def curie_target(token):
    """'natTi' -> 'natTI', '27Al' -> '27AL', '60Ni' -> '60NI'."""
    token = token.strip()
    if token.lower().startswith('nat'):
        return 'nat' + token[3:].upper()
    return curie_nuclide(token)
