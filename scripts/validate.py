"""Build-time validation for the generated databases.

Each database gets a list of named checks; any failure blocks release.
Checks compare against the schema contract curie's readers require and,
where a previous-generation file is available (CURIE_V1_DIR, default
~/.local/share/curie), produce an old-vs-new coverage diff as the review
exhibit.

Usage: python3 validate.py build/iaea_monitors.db
"""
import os
import sqlite3
import sys
from pathlib import Path

V1_DIR = Path(os.environ.get('CURIE_V1_DIR', Path.home() / '.local/share/curie'))

PASS, FAIL = 'pass', 'FAIL'
_results = []


def check(name, ok, detail=''):
    _results.append((PASS if ok else FAIL, name, detail))
    print(f'  [{PASS if ok else FAIL}] {name}' + (f' — {detail}' if detail else ''))
    return ok


def q(con, sql):
    return con.execute(sql).fetchall()


# --- iaea_monitors ---

IAEA_PARTICLES = {'p': 'proton', 'd': 'deuteron', 'h': 'helion', 'a': 'alpha'}
IAEA_DATA_COLS = ['target', 'outgoing', 'product', 'energy',
                  'cross_section', 'unc_cross_section']


def validate_iaea_monitors(con):
    tables = {r[0] for r in q(con, "SELECT name FROM sqlite_master WHERE type='table'")}
    check('tables', tables >= set(IAEA_PARTICLES.values()) | {'all_reactions', '_version'},
          str(sorted(tables)))
    check('legacy neutron/gamma tables absent', not ({'neutron', 'gamma'} & tables))

    cols = [c[1] for c in q(con, 'PRAGMA table_info(all_reactions)')]
    check('all_reactions column order (product last = retrieval label)',
          cols == ['target', 'incident', 'outgoing', 'product'], str(cols))
    for t in IAEA_PARTICLES.values():
        cols = [c[1] for c in q(con, f'PRAGMA table_info({t})')]
        check(f'{t} columns', cols == IAEA_DATA_COLS, str(cols))

    # uniqueness (Reaction._check requires exactly one match)
    dup = q(con, 'SELECT target,incident,outgoing,product,COUNT(*) c FROM all_reactions '
                 'GROUP BY 1,2,3,4 HAVING c>1')
    check('no duplicate index rows', not dup, str(dup[:3]))

    # bidirectional index <-> data joins
    missing_data, missing_index = [], []
    for inc, t in IAEA_PARTICLES.items():
        idx = set(q(con, f"SELECT target,outgoing,product FROM all_reactions WHERE incident='{inc}'"))
        dat = set(q(con, f'SELECT DISTINCT target,outgoing,product FROM {t}'))
        missing_data += [(t, r) for r in idx - dat]
        missing_index += [(t, r) for r in dat - idx]
    check('every index row has data', not missing_data, str(missing_data[:3]))
    check('every data reaction is indexed', not missing_index, str(missing_index[:3]))

    # grids: strictly sorted, no duplicate abscissae per reaction
    bad_sort = []
    for inc, t in IAEA_PARTICLES.items():
        for tgt, out, prod in q(con, f'SELECT DISTINCT target,outgoing,product FROM {t}'):
            e = [r[0] for r in q(con, f"SELECT energy FROM {t} WHERE target='{tgt}' "
                                      f"AND outgoing='{out}' AND product='{prod}' ORDER BY rowid")]
            if any(b <= a for a, b in zip(e, e[1:])):
                bad_sort.append((t, tgt, out, prod))
    check('grids strictly increasing (sorted, no duplicate abscissae)',
          not bad_sort, str(bad_sort[:3]))

    # grid coverage: the stored grids must equal a fresh re-parse of the
    # mirrored sources exactly — catches any silent row loss or truncation
    # (e.g. rows whose uncertainty cell is blank must still be present)
    sys.path.insert(0, str(Path(__file__).parent / 'build'))
    import iaea_monitors as builder
    reactions, _, _ = builder.collect()
    mismatch = []
    for key, info in reactions.items():
        tgt, inc, out, prod = key
        want = builder.build_grid(key, info, announce=None)
        got = q(con, f"SELECT energy,cross_section,unc_cross_section FROM "
                     f"{IAEA_PARTICLES[inc]} WHERE target='{tgt}' AND "
                     f"outgoing='{out}' AND product='{prod}' ORDER BY energy")
        if [tuple(r) for r in got] != want:
            detail = (f'{len(got)} rows to {got[-1][0] if got else 0} MeV vs '
                      f'source {len(want)} rows to {want[-1][0]} MeV')
            mismatch.append(('%s(%s,%s)%s' % key, detail))
    check(f'stored grids identical to source re-parse ({len(reactions)} reactions)',
          not mismatch, str(mismatch[:3]))

    # physical sanity
    neg = sum(q(con, f'SELECT COUNT(*) FROM {t} WHERE cross_section<0 OR unc_cross_section<0')[0][0]
              for t in IAEA_PARTICLES.values())
    check('no negative values', neg == 0, f'{neg} rows')

    # spot anchors (known defects of the previous generation must be gone)
    spike = q(con, "SELECT energy,cross_section FROM proton WHERE target='natFE' "
                   "AND product='51CR' AND cross_section>1000")
    check('natFE(p,x)51CR spike absent (no sigma > 1000 mb)', not spike, str(spike))
    ti = q(con, "SELECT COUNT(*),MAX(cross_section) FROM proton WHERE target='natTI' AND product='48V'")
    check('natTi(p,x)48V present, peak sigma plausible (300-500 mb)',
          ti[0][0] > 0 and 300 < ti[0][1] < 500, str(ti[0]))

    ver = q(con, 'SELECT library,generation,build_date FROM _version')
    check('_version stamp', len(ver) == 1 and ver[0][0] == 'iaea_monitors'
          and ver[0][1] == 'v2', str(ver))


def coverage_diff_iaea(con):
    v1p = V1_DIR / 'iaea_monitors.db'
    if not v1p.exists():
        print('  (no v1 file for coverage diff)')
        return
    v1 = sqlite3.connect(v1p)
    old = set(q(v1, 'SELECT target,incident,outgoing,product FROM all_reactions'))
    new = set(q(con, 'SELECT target,incident,outgoing,product FROM all_reactions'))
    added, dropped = sorted(new - old), sorted(old - new)
    print(f'\n  coverage vs v1: {len(old)} -> {len(new)} reactions; '
          f'+{len(added)} added, -{len(dropped)} dropped')
    for r in added:
        print('    + %s(%s,%s)%s' % r)
    for r in dropped:
        print('    - %s(%s,%s)%s' % r)
    # energy-range regressions on reactions present in both generations
    # (report-only: a narrowed evaluation is legitimate but must be seen)
    for tgt, inc, out, prod in sorted(old & new):
        t = IAEA_PARTICLES.get(inc)
        if t is None:
            continue
        args = (tgt, out, prod)
        sql = (f"SELECT MAX(energy) FROM {t} WHERE target=? AND outgoing=? "
               f"AND product=?")
        e1 = v1.execute(sql, args).fetchone()[0]
        e2 = con.execute(sql, args).fetchone()[0]
        if e1 and e2 and e2 < e1 - 1e-9:
            print(f'    range regression: {tgt}({inc},{out}){prod} '
                  f'max E {e1} -> {e2} MeV')


VALIDATORS = {
    'iaea_monitors.db': (validate_iaea_monitors, coverage_diff_iaea),
}


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(f'usage: validate.py <db>  (one of: {", ".join(VALIDATORS)})')
    path = Path(sys.argv[1])
    if not path.exists():
        sys.exit(f'no such file: {path}')
    if path.name not in VALIDATORS:
        sys.exit(f'no validator for {path.name} (one of: {", ".join(VALIDATORS)})')
    con = sqlite3.connect(path)
    validate, diff = VALIDATORS[path.name]
    print(f'validating {path}')
    validate(con)
    diff(con)
    n_fail = sum(1 for r in _results if r[0] == FAIL)
    print(f'\n{len(_results)} checks, {n_fail} failures')
    sys.exit(1 if n_fail else 0)
