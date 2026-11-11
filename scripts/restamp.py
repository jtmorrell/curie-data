"""Set the generation stamp of a published database without rebuilding it.

    python scripts/restamp.py <source.db> <out.db> <generation> --sha256 <hex>

Copies the source byte for byte, then sets the one `generation` cell of its
`_version` table. No other cell changes and no VACUUM runs, so the result
differs from the source only in that cell (and the page and header bytes
SQLite rewrites to store it), and every other table stays equal by content.

The source's sha256 must be given and must match, so the stamp is only ever
applied to the exact file that was published. The result's sha256 and size
are printed for the registry and the data manifest.
"""
import argparse
import hashlib
import os
import shutil
import sqlite3
import sys
from pathlib import Path


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _version_rows(path):
    """The `_version` rows of the database at path, opened read-only."""
    con = sqlite3.connect(f'file:{Path(path).resolve()}?mode=ro', uri=True)
    try:
        if not con.execute("SELECT 1 FROM sqlite_master WHERE type='table' "
                           "AND name='_version'").fetchone():
            return None
        return con.execute('SELECT * FROM _version').fetchall()
    finally:
        con.close()


def restamp(source, out, generation, sha256):
    """Write a copy of source to out with `_version.generation` set to
    generation, and return (sha256, size) of out.

    Raises ValueError, writing nothing, when the source's sha256 is not
    sha256, when out is the source, or when the source has no one-row
    `_version` table.
    """
    source, out = Path(source), Path(out)
    actual = sha256_of(source)
    if actual != sha256.strip().lower():
        raise ValueError(f'{source}: sha256 is {actual}, not the given '
                         f'{sha256}')
    if out.resolve() == source.resolve():
        raise ValueError(f'{out}: the output must not be the source')
    rows = _version_rows(source)
    if rows is None:
        raise ValueError(f'{source}: no _version table')
    if len(rows) != 1:
        raise ValueError(f'{source}: _version holds {len(rows)} rows, not 1')

    tmp = out.with_name(out.name + '.restamp-tmp')
    try:
        shutil.copyfile(source, tmp)
        con = sqlite3.connect(tmp)
        try:
            con.execute('UPDATE _version SET generation = ?', (generation,))
            con.commit()
        finally:
            con.close()
        os.replace(tmp, out)
    finally:
        if tmp.exists():
            tmp.unlink()
    return sha256_of(out), out.stat().st_size


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('source', help='the published database to restamp')
    p.add_argument('out', help='where to write the restamped copy')
    p.add_argument('generation', help='the new _version.generation, e.g. v2.1')
    p.add_argument('--sha256', required=True,
                   help="the source's expected sha256 (refuses on mismatch)")
    args = p.parse_args(argv)
    try:
        digest, size = restamp(args.source, args.out, args.generation,
                               args.sha256)
    except ValueError as e:
        print(f'refused: {e}', file=sys.stderr)
        return 1
    print(f'sha256 {digest}')
    print(f'size {size}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
