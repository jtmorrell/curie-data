"""Tests for scripts/restamp.py and the decay build's generation stamp.

Each test builds a small SQLite file shaped like a built database (a one-row
`_version` table beside ordinary tables) and restamps it; nothing reads the
shipped data.

Run with `python3 -m pytest tests/`.
"""
import ast
import hashlib
import inspect
import sqlite3
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'scripts' / 'build'))
import restamp  # noqa: E402

VERSION_ROW = ('decay', 'v2', 'a source', 'source 1.0', '2026-07-16', 'abc1234')


def _make_db(path, with_version=True):
    con = sqlite3.connect(path)
    con.execute('CREATE TABLE chart (isotope TEXT, half_life REAL, '
                'unc_half_life REAL)')
    con.executemany('INSERT INTO chart VALUES (?,?,?)',
                    [('152EUg', 4.2716e8, 1.6e5), ('60COg', 1.6634e8, 8.6e3),
                     ('137CSg', 9.4925e8, None)])
    con.execute('CREATE TABLE coincidences (isotope TEXT, e1 REAL, e2 REAL, '
                'intensity REAL)')
    con.executemany('INSERT INTO coincidences VALUES (?,?,?,?)',
                    [('60COg', 1173.228, 1332.492, 0.9985),
                     ('152EUg', 121.7817, 344.2785, 0.0731)])
    con.execute('CREATE INDEX coinc_iso ON coincidences (isotope)')
    if with_version:
        con.execute('CREATE TABLE _version (library TEXT, generation TEXT, '
                    'source TEXT, source_version TEXT, build_date TEXT, '
                    'pipeline TEXT)')
        con.execute('INSERT INTO _version VALUES (?,?,?,?,?,?)', VERSION_ROW)
    con.commit()
    con.close()
    return path


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _contents(path):
    """Schema and every table's rows in rowid order."""
    con = sqlite3.connect(path)
    schema = con.execute('SELECT * FROM sqlite_master ORDER BY rowid').fetchall()
    tables = {name: con.execute(f'SELECT * FROM "{name}" ORDER BY rowid').fetchall()
              for (name,) in con.execute(
                  "SELECT name FROM sqlite_master WHERE type='table'")}
    con.close()
    return schema, tables


def test_restamp_changes_only_the_generation(tmp_path, capsys):
    src = _make_db(tmp_path / 'src.db')
    src_sha = _sha(src)
    out = tmp_path / 'out.db'
    assert restamp.main([str(src), str(out), 'v2.1', '--sha256', src_sha]) == 0

    schema0, tables0 = _contents(src)
    schema1, tables1 = _contents(out)
    assert tables1['_version'] == [('decay', 'v2.1') + VERSION_ROW[2:]]
    assert tables0['_version'] == [VERSION_ROW]
    assert schema1 == schema0
    assert {k: v for k, v in tables1.items() if k != '_version'} == \
           {k: v for k, v in tables0.items() if k != '_version'}
    assert out.stat().st_size == src.stat().st_size
    assert _sha(src) == src_sha                      # source untouched
    printed = capsys.readouterr().out.split('\n')
    assert f'sha256 {_sha(out)}' in printed
    assert f'size {out.stat().st_size}' in printed
    assert not list(tmp_path.glob('*restamp-tmp*'))


def test_wrong_sha256_refuses_and_writes_nothing(tmp_path):
    src = _make_db(tmp_path / 'src.db')
    out = tmp_path / 'out.db'
    assert restamp.main([str(src), str(out), 'v2.1', '--sha256', '0' * 64]) == 1
    assert sorted(p.name for p in tmp_path.iterdir()) == ['src.db']
    with pytest.raises(ValueError, match='sha256'):
        restamp.restamp(src, out, 'v2.1', '0' * 64)
    assert not out.exists()


def test_sha256_is_required(tmp_path):
    src = _make_db(tmp_path / 'src.db')
    with pytest.raises(SystemExit):
        restamp.main([str(src), str(tmp_path / 'out.db'), 'v2.1'])
    assert sorted(p.name for p in tmp_path.iterdir()) == ['src.db']


def test_source_without_version_refuses_and_writes_nothing(tmp_path):
    src = _make_db(tmp_path / 'src.db', with_version=False)
    out = tmp_path / 'out.db'
    assert restamp.main([str(src), str(out), 'v2.1',
                         '--sha256', _sha(src)]) == 1
    assert sorted(p.name for p in tmp_path.iterdir()) == ['src.db']


def test_output_equal_to_source_refuses(tmp_path):
    src = _make_db(tmp_path / 'src.db')
    sha = _sha(src)
    assert restamp.main([str(src), str(src), 'v2.1', '--sha256', sha]) == 1
    assert _sha(src) == sha


def test_decay_build_writes_its_own_generation():
    import common
    import decay
    assert decay.DECAY_GENERATION == 'v2.1'
    assert inspect.signature(common.write_version).parameters[
        'generation'].default == 'v2'
    # the build's write_version call passes the constant
    tree = ast.parse((ROOT / 'scripts' / 'build' / 'decay.py').read_text())
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
             and getattr(n.func, 'id', None) == 'write_version']
    assert len(calls) == 1
    kw = {k.arg: k.value for k in calls[0].keywords}
    assert isinstance(kw.get('generation'), ast.Name)
    assert kw['generation'].id == 'DECAY_GENERATION'
