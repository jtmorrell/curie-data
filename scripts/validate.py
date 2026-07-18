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
import zipfile
from collections import defaultdict
from pathlib import Path

V1_DIR = Path(os.environ.get('CURIE_V1_DIR', Path.home() / '.local/share/curie'))

PASS, FAIL = 'pass', 'FAIL'
_results = []


def _chart_names():
    """Isotope names from the decay database, for product-resolution checks."""
    p = V1_DIR / 'decay.db'
    if not p.exists():
        return None
    con = sqlite3.connect(p)
    return {r[0] for r in con.execute('SELECT name FROM chart')}


def resolvable_products(con, particle_letters=False, blocking=True):
    """Every all_reactions.product must resolve to a decay-database isotope
    (after curie's normalization: bare names are ground states, bare 'm' is
    the first isomer). Empty products and — where allowed — particle
    letters pass. This is the gate that catches corrupted or mis-numbered
    product names at build time."""
    chart = _chart_names()
    if chart is None:
        print('  (no decay.db for product-resolution check)')
        return
    import re as _re
    bad = []
    for (p,) in q(con, 'SELECT DISTINCT product FROM all_reactions'):
        # the particle letters include 'n' and 'g' (the nprod/gprod
        # production-sum specials store their emitted particle as product)
        if not p or (particle_letters and p in ('n', 'g', 'p', 'd', 't', 'h', 'a')):
            continue
        # curie's normalization: bare 'm' means the first isomer; a name
        # with no state suffix means the ground state
        if p.endswith('m'):
            name = p + '1'
        elif _re.search(r'(m\d+|g)$', p):
            name = p
        else:
            name = p + 'g'
        if name not in chart:
            bad.append(p)
    label = 'products resolve in decay.db'
    if blocking:
        check(label, not bad, str(bad[:5]))
    else:
        print(f'  [info] {label}: {len(bad)} unresolved '
              f'(report-only until the rebuilt decay database ships): {bad[:8]}')


def check(name, ok, detail=''):
    _results.append((PASS if ok else FAIL, name, detail))
    print(f'  [{PASS if ok else FAIL}] {name}' + (f' — {detail}' if detail else ''))
    return ok


def q(con, sql, args=()):
    return con.execute(sql, args).fetchall()


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
    # collect() discovers reactions by glob, which returns nothing on a
    # missing or mispointed archive — an empty re-parse must fail loudly,
    # not let the identity check below pass vacuously
    check('source archive present (re-parse found reactions)',
          len(reactions) > 0, str(builder.ARCHIVE_ROOT / 'iaea-medical'))
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


# --- IRDFF ---

def validate_irdff(con):
    idx = q(con, 'SELECT target,outgoing,product,table_name FROM all_reactions')
    check('reaction count', len(idx) == 119, str(len(idx)))
    cols = [c[1] for c in q(con, 'PRAGMA table_info(all_reactions)')]
    check('all_reactions column order (table_name last = retrieval label)',
          cols == ['target', 'outgoing', 'product', 'table_name'], str(cols))
    bad_name = [(t, o, p, n) for t, o, p, n in idx if n != f'_{t}_n_{o}_{p}']
    check('table names derive from reaction identity', not bad_name, str(bad_name[:3]))

    tables = {r[0] for r in q(con, "SELECT name FROM sqlite_master WHERE type='table'")}
    named = {n for *_, n in idx}
    check('index <-> tables bidirectional',
          named <= tables and tables - named == {'all_reactions', '_version'},
          str(sorted(tables - named - {'all_reactions', '_version'})[:3]))

    bad_cols, bad_sort, empty = [], [], []
    for t in sorted(named):
        cols = [c[1] for c in q(con, f'PRAGMA table_info({t})')]
        if cols != ['energy', 'cross_section', 'unc_cross_section']:
            bad_cols.append(t)  # positional SELECT * on read: order load-bearing
        e = [r[0] for r in q(con, f'SELECT energy FROM {t} ORDER BY rowid')]
        if not e:
            empty.append(t)
        if any(b < a for a, b in zip(e, e[1:])):
            bad_sort.append(t)  # nondecreasing; duplicate-E step points are meaningful
    check('data-table column order (positional read)', not bad_cols, str(bad_cols[:3]))
    check('grids nondecreasing and non-empty', not bad_sort and not empty,
          str((bad_sort + empty)[:3]))

    th = q(con, 'SELECT cross_section FROM _115IN_n_g_116INm '
                'WHERE energy BETWEEN 0.02529 AND 0.02531')
    check('115IN(n,g)116INm thermal anchor ~162 b',
          len(th) == 1 and 150 < th[0][0] < 175, str(th))

    # blocking identity gate: stored data must equal a fresh source re-parse
    sys.path.insert(0, str(Path(__file__).parent / 'build'))
    import irdff as builder
    blocks = builder.parse_tab(
        zipfile.ZipFile(builder.ARCHIVE_ROOT / 'irdff-ii' / 'IRDFF-II_TAB.zip')
        .read('IRDFF-II.tab.txt').decode('utf-8', 'replace'))
    mismatch = []
    for t, o, p, n in idx:
        src = blocks.get((t, o, p)) or blocks.get((t, o, ''))
        got = [tuple(r) for r in q(con, f'SELECT * FROM {n} ORDER BY rowid')]
        if src is None or got != src:
            mismatch.append(n)
    check('stored data identical to source re-parse (119 reactions)',
          not mismatch, str(mismatch[:3]))

    resolvable_products(con)
    ver = q(con, 'SELECT library,generation FROM _version')
    check('_version stamp', ver == [('IRDFF', 'v2')], str(ver))


def coverage_diff_irdff(con):
    v1p = V1_DIR / 'IRDFF.db'
    if not v1p.exists():
        print('  (no v1 file for coverage diff)')
        return
    v1 = sqlite3.connect(v1p)
    sys.path.insert(0, str(Path(__file__).parent / 'build'))
    from irdff import V1_NAME_FIXES
    ident, renamed, differ = 0, [], []
    v2names = {(t, o, p): n for t, o, p, n in
               q(con, 'SELECT target,outgoing,product,table_name FROM all_reactions')}
    for t, o, p, n in q(v1, 'SELECT target,outgoing,product,table_name FROM all_reactions'):
        n2 = v2names.get((t, o, V1_NAME_FIXES.get(p, p)))
        if n2 is None:
            differ.append((n, 'MISSING'))
            continue
        r1 = q(v1, f'SELECT * FROM {n} ORDER BY energy')
        r2 = q(con, f'SELECT * FROM {n2} ORDER BY energy')
        if r1 == r2:
            ident += 1
            if n != n2:
                renamed.append(f'{n} -> {n2}')
        else:
            differ.append((n, f'{len(r1)} vs {len(r2)} rows'))
    print(f'\n  vs v1: {ident}/119 data-identical; renamed: {renamed}; '
          f'differing: {differ[:5] if differ else "none"}')


# --- endf ---

def _endf_table(target):
    import re as _re
    if target.startswith('nat'):
        # natural-element targets: 'natFE' -> FE_000 (elemental A=000);
        # the generic rule below would strip to 'FE_' and misread the
        # 'm' of a natSM as an isomer flag
        return target[3:] + '_000'
    return (''.join(_re.findall('[A-Z]+', target)) + '_' +
            ''.join(_re.findall('[0-9]+', target)) +
            ('m' if 'm' in target else ''))


def validate_endf(con):
    idx = q(con, 'SELECT target,outgoing,product,label FROM all_reactions')
    tables = {r[0] for r in q(con, "SELECT name FROM sqlite_master WHERE type='table'")}
    missing_tab = [t for t in sorted({r[0] for r in idx})
                   if _endf_table(t) not in tables]
    check('every target table derivable and present (shard round-trip)',
          not missing_tab, str(missing_tab[:3]))
    targets = _wide_table_checks(con, idx,
                                 ['target', 'outgoing', 'product', 'label'])
    check('material count (VIII.1 = 558 minus the excluded free-neutron '
          'evaluation)', len(targets) == 557, str(len(targets)))

    # spot anchors at 0.0253 eV, incl. the computed-absorption column and
    # residual-product identities (the two derived quantities)
    for tab, col, lo, hi in [('IN_115', '_g_', 195, 210),
                             ('FE_56', '_g_', 2.5, 2.65),
                             ('U_235', '_fis_', 570, 600),
                             ('AU_197', '_g_', 96, 101),
                             ('B_10', '_abs_', 3800, 3900),
                             ('U_235', '_abs_', 95, 105)]:
        r = q(con, f'SELECT energy,"{col}" FROM {tab} WHERE '
                   f'energy BETWEEN 0.02 AND 0.03 ORDER BY ABS(energy-0.0253) LIMIT 1')
        ok = bool(r) and lo < r[0][1] < hi
        check(f'{tab} {col} thermal anchor in ({lo}, {hi}) b', ok, str(r))
    prod_anchors = [('56FE', 'p', '56MN'), ('27AL', 'a', '24NA'),
                    ('115IN', 'g', '116IN'), ('139LA', '5n', '135LA')]
    bad_prod = [(t, o, want) for t, o, want in prod_anchors
                if q(con, f"SELECT product FROM all_reactions WHERE "
                          f"target='{t}' AND outgoing='{o}'") != [(want,)]]
    check('residual-product identity anchors', not bad_prod, str(bad_prod))

    resolvable_products(con, particle_letters=True, blocking=False)
    ver = q(con, 'SELECT library,generation FROM _version')
    check('_version stamp', ver == [('endf', 'v2')], str(ver))


def coverage_diff_endf(con):
    v1p = V1_DIR / 'endf.db'
    if not v1p.exists():
        print('  (no v1 file for coverage diff)')
        return
    v1 = sqlite3.connect(v1p)
    old = {(t, o) for t, o in q(v1, 'SELECT target,outgoing FROM all_reactions')}
    new = {(t, o) for t, o in q(con, 'SELECT target,outgoing FROM all_reactions')}
    old_t = {t for t, _ in old}
    new_t = {t for t, _ in new}
    added_t, dropped_t = sorted(new_t - old_t), sorted(old_t - new_t)
    print(f'\n  materials: {len(old_t)} -> {len(new_t)} '
          f'(+{len(added_t)} added, -{len(dropped_t)} dropped)')
    print(f'  new materials: {added_t[:12]}{"..." if len(added_t) > 12 else ""}')
    if dropped_t:
        print(f'  DROPPED materials: {dropped_t}')
    both = old_t & new_t
    ch_lost = sorted((t, o) for t, o in old if t in both and (t, o) not in new)
    ch_new = [(t, o) for t, o in new if t in both and (t, o) not in old]
    print(f'  channels on common materials: +{len(ch_new)} added, '
          f'-{len(ch_lost)} dropped')
    for t, o in ch_lost[:20]:
        print(f'    - {t}(n,{o})')
    # per-channel energy-range regressions on common (target, outgoing)
    # pairs (report-only: VIII.1 grids legitimately differ from VII.1, but
    # a shortened range must be seen)
    regress = []
    for t, o in sorted(old & new):
        lab = q(con, f"SELECT label FROM all_reactions WHERE target='{t}' AND outgoing='{o}'")
        lab1 = q(v1, f"SELECT label FROM all_reactions WHERE target='{t}' AND outgoing='{o}'")
        if len(lab) != 1 or len(lab1) != 1:
            continue
        tab = _endf_table(t)
        e2 = q(con, f'SELECT MAX(energy) FROM {tab} WHERE "{lab[0][0]}" > 0')[0][0]
        e1 = q(v1, f'SELECT MAX(energy) FROM {tab} WHERE "{lab1[0][0]}" > 0')[0][0]
        if e1 and e2 and e2 < 0.5 * e1:
            regress.append((t, o, e1 / 1e6, e2 / 1e6))
    print(f'  channel range regressions >2x on common pairs: {len(regress)}')
    for t, o, a, b in regress[:15]:
        print(f'    {t}(n,{o}): {a:.3g} -> {b:.3g} MeV')


# --- tendl (exclusive + residual-production) ---

def _wide_table_checks(con, idx, index_cols, allow_negative=()):
    """Shared structural and value checks for the wide-table libraries
    (endf and the five tendl databases): index/table label joins, strictly
    increasing non-empty grids, and a full value sweep — no NULL and no
    negative entry in the energy column or any channel column, except
    columns named in `allow_negative` (documented source artifacts)."""
    cols = [c[1] for c in q(con, 'PRAGMA table_info(all_reactions)')]
    check('all_reactions column order (label last)', cols == index_cols, str(cols))
    key = ','.join(index_cols[:-1])
    dup = q(con, f'SELECT {key},COUNT(*) c FROM all_reactions '
                 f'GROUP BY {key} HAVING c>1')
    check('no ambiguous index rows', not dup, str(dup[:3]))
    targets = sorted({r[0] for r in idx})
    labels_by_target = defaultdict(set)
    for r in idx:
        labels_by_target[r[0]].add(r[-1])
    bad_join, bad_sort, bad_val = [], [], []
    for tgt in targets:
        tab = _endf_table(tgt)
        have = {c[1] for c in q(con, f'PRAGMA table_info({tab})')} - {'energy'}
        if labels_by_target[tgt] != have:
            bad_join.append((tab, sorted(labels_by_target[tgt] ^ have)[:4]))
        e = [r[0] for r in q(con, f'SELECT energy FROM {tab} ORDER BY rowid')]
        if not e or any(b <= a for a, b in zip(e, e[1:])):
            bad_sort.append(tab)
        # value sweep, one scan per <=400-column chunk: MIN and NULL count
        # for the energy column and every channel column
        sweep_cols = ['energy'] + sorted(have)
        for i in range(0, len(sweep_cols), 400):
            chunk = sweep_cols[i:i + 400]
            expr = ', '.join(f'MIN("{c}"), COUNT(*)-COUNT("{c}")' for c in chunk)
            row = q(con, f'SELECT {expr} FROM {tab}')[0]
            for c, mn, n_null in zip(chunk, row[0::2], row[1::2]):
                if mn is None or n_null or (mn < 0 and c not in allow_negative):
                    bad_val.append((tab, c, (mn, n_null)))
    check('index labels == table columns (bidirectional)', not bad_join, str(bad_join[:3]))
    check('grids strictly increasing and non-empty', not bad_sort, str(bad_sort[:3]))
    check('no NULL or negative values in energy or any channel', not bad_val,
          str(bad_val[:3]))
    return targets


def validate_tendl(con):
    idx = q(con, 'SELECT target,outgoing,product,label FROM all_reactions')
    # _el_ may go negative in isolated rows at the bottom of the resonance
    # region: the ~500-point TALYS grid interpolates through resonance dips
    # and the source prints the artifact (3 targets in TENDL-2025; the
    # previous generation shipped one as well). Stored source-faithfully —
    # the resonance region is documented as non-resolving and 'best' serves
    # it from irdff/endf.
    targets = _wide_table_checks(con, idx,
                                 ['target', 'outgoing', 'product', 'label'],
                                 allow_negative=('_el_',))
    check('target count (curated cut incl. isomeric targets, neutron)',
          350 < len(targets) < 450, str(len(targets)))

    # source-identity gate: stored content equals a fresh table re-parse
    # (the mirror may hold targets outside the curated cut — skip those)
    sys.path.insert(0, str(Path(__file__).parent / 'build'))
    import tendl_neutron as builder
    src = builder.ARCHIVE_ROOT / 'tendl-2025' / 'neutron_files'
    cut = builder.curated_targets()
    mismatch = []
    for tdir in sorted(t for el in src.iterdir() if el.is_dir()
                       for t in el.iterdir() if (t / 'tables' / 'xs').is_dir()
                       and builder.mirror_target(t.name) in cut):
        tgt, tab, cols, grid, data = builder.assemble(tdir)
        if tab is None:
            continue
        got = q(con, f'SELECT * FROM {tab} ORDER BY rowid')
        if len(got) != len(data) or [tuple(r) for r in got] != [tuple(r) for r in data.tolist()]:
            mismatch.append(tab)
    check('stored tables identical to source re-parse', not mismatch, str(mismatch[:3]))

    # NOTE: no thermal anchor — the TALYS tables' ~500-point grid does not
    # resolve the resonance region (the previous generation reads 40 b at
    # the 115In thermal point against the true 202 b), and 'best' serves
    # thermal neutrons from irdff/endf. Anchor on the fast region instead.
    r = q(con, 'SELECT energy,"_2n_" FROM IN_115 ORDER BY ABS(energy-14.0) LIMIT 1')
    check('IN_115 (n,2n) 14 MeV anchor in (1000, 1900) mb',
          bool(r) and 1000 < r[0][1] < 1900, str(r))
    resolvable_products(con, particle_letters=True, blocking=False)
    ver = q(con, 'SELECT library,generation FROM _version')
    check('_version stamp', ver == [('tendl', 'v2')], str(ver))


def coverage_diff_tendl(con):
    v1p = V1_DIR / 'tendl.db'
    if not v1p.exists():
        print('  (no v1 file for coverage diff)')
        return
    v1 = sqlite3.connect(v1p)
    old = {(t, o) for t, o in q(v1, 'SELECT target,outgoing FROM all_reactions')}
    new = {(t, o) for t, o in q(con, 'SELECT target,outgoing FROM all_reactions')}
    old_t, new_t = {t for t, _ in old}, {t for t, _ in new}
    print(f'\n  targets: {len(old_t)} -> {len(new_t)} (+{len(new_t - old_t)}, '
          f'-{len(old_t - new_t)}); dropped: {sorted(old_t - new_t) or "none"}')
    both = old_t & new_t
    lost = [(t, o) for t, o in old if t in both and (t, o) not in new]
    lvl = [x for x in lost if x[1] and x[1][-1].isdigit()]
    print(f'  channels on common targets: +{len([x for x in new if x[0] in both and x not in old])} '
          f'added, -{len(lost)} dropped ({len(lvl)} discrete-level entries no '
          f'longer published; {len(lost) - len(lvl)} other)')
    other = [x for x in lost if x not in lvl]
    for t, o in other[:10]:
        print(f'    - {t}(n,{o})')


def _validate_tendl_rp(con, beam, letter):
    idx = q(con, 'SELECT target,product,label FROM all_reactions')
    targets = _wide_table_checks(con, idx, ['target', 'product', 'label'])
    check(f'target count (curated cut + natural-element targets, {beam})',
          400 < len(targets) < 560, str(len(targets)))
    bad_suffix = sorted({p for _, p, _ in idx
                         if not (p.endswith('g') or 'm' in p[-3:])})
    check('every product carries an explicit state suffix', not bad_suffix,
          str(bad_suffix[:5]))

    # source-identity gate: stored content equals the builder's own
    # assembly of a fresh source re-parse (the same code path the build ran)
    sys.path.insert(0, str(Path(__file__).parent / 'build'))
    import tendl_rp as builder
    src = builder.ARCHIVE_ROOT / 'tendl-2025' / f'{beam}_files'
    # the mirror may hold targets outside the curated cut — skip those
    cut = builder.curated_targets()
    mismatch = []
    for tdir in sorted(t for el in src.iterdir() if el.is_dir()
                       for t in el.iterdir()
                       if (t / 'tables' / 'residual').is_dir()
                       and builder.mirror_target(t.name) in cut):
        tgt, tab, cols, grid, data = builder.assemble(tdir)
        if tab is None:
            continue
        got = q(con, f'SELECT * FROM {tab} ORDER BY rowid')
        if [tuple(r) for r in got] != [tuple(r) for r in data.tolist()]:
            mismatch.append(tab)
    check('stored tables identical to source re-parse', not mismatch, str(mismatch[:3]))

    # natural-element tables (abundance-weighted isotopic sums): present
    # exactly for the fully covered elements, and byte-identical to a
    # fresh recomputation through the builder's shared recipe
    tables = {r[0] for r in q(con, "SELECT name FROM sqlite_master WHERE type='table'")}
    nat_wrong, nat_bad = [], []
    for el, parts in sorted(builder.nat_constituents().items()):
        complete = all(tab in tables for _, tab, _ in parts)
        if (f'{el}_000' in tables) != complete:
            nat_wrong.append(el)
            continue
        if not complete:
            continue
        cols, grid, data = builder.assemble_nat(con, parts)
        names = [c[1] for c in q(con, f'PRAGMA table_info({el}_000)')][1:]
        got = q(con, f'SELECT * FROM {el}_000 ORDER BY rowid')
        if names != [f'_{c}' for c in cols] or \
                [tuple(r) for r in got] != [tuple(r) for r in data.tolist()]:
            nat_bad.append(el)
    check('nat table present iff element fully covered', not nat_wrong,
          str(nat_wrong[:5]))
    check('nat tables identical to weighted-sum recomputation', not nat_bad,
          str(nat_bad[:3]))

    if beam == 'proton':
        _tendl_rp_endf_checks(con, builder, src)
        _nat_monitor_plausibility(con)

    resolvable_products(con, blocking=False)
    ver = q(con, 'SELECT library,generation FROM _version')
    check('_version stamp', ver == [(f'tendl_{letter}_rp', 'v2')], str(ver))


def _nat_monitor_plausibility(con):
    """natTi(p,x)48V against the IAEA monitor evaluation — a plausibility
    band, not a truth anchor: TENDL is TALYS systematics while the IAEA
    monitors are fit to data, so only peak position (3 MeV) and magnitude
    (factor of 2) are gated; the ratio is reported for the record."""
    import numpy as np
    iaea_path = Path(__file__).resolve().parents[1] / 'build' / 'iaea_monitors.db'
    if not iaea_path.exists():
        check('nat-vs-monitor plausibility (iaea_monitors.db present)',
              False, str(iaea_path))
        return
    icon = sqlite3.connect(iaea_path)
    ref = np.array(icon.execute(
        "SELECT energy, cross_section FROM proton WHERE target='natTI' "
        "AND product='48V' ORDER BY energy").fetchall())
    icon.close()
    got = np.array(q(con, 'SELECT energy, "_48Vg" FROM TI_000 ORDER BY rowid'))
    ri, gi = ref[:, 1].argmax(), got[:, 1].argmax()
    de = abs(got[gi, 0] - ref[ri, 0])
    ratio = got[gi, 1] / ref[ri, 1]
    print(f'  natTI(p,x)48V vs IAEA monitor: peak {got[gi, 1]:.0f} mb @ '
          f'{got[gi, 0]:.1f} MeV vs {ref[ri, 1]:.0f} mb @ {ref[ri, 0]:.1f} '
          f'MeV (ratio {ratio:.2f})')
    check('natTI(p,x)48V within the monitor plausibility band '
          '(peak +-3 MeV, factor 2)', de <= 3.0 and 0.5 <= ratio <= 2.0,
          f'dE={de:.1f} MeV, ratio={ratio:.2f}')


def _tendl_rp_endf_checks(con, builder, src):
    """The proton database carries ENDF-6-ingested targets (the TENDL-2023
    Nb exception). Two gates: the extractor must reproduce the TENDL-2025
    Fe056 tables from the same generation's s30 file (proving the
    extraction physics), and each ingested table must equal a fresh
    re-extraction (the same identity standard as the mirror targets)."""
    import re as _re
    import numpy as np
    from tendl_rp_endf import assemble_endf, endf_targets, extract_products

    fe_file = builder.ARCHIVE_ROOT / 'tendl-2025' / 'p-Fe056-s30.tendl'
    fe_tables = builder.ARCHIVE_ROOT / 'tendl-2025' / 'proton_files' / 'Fe' / 'Fe056' / 'tables' / 'residual'
    if not (fe_file.exists() and fe_tables.is_dir()):
        check('rp-from-ENDF-6 identity material present (Fe056 s30 + tables)',
              False, str(fe_file))
        return
    ext = extract_products(fe_file)
    ext_nuclides = {_re.sub(r'(g|m\d+)$', '', k) for k in ext}
    bad, compared, state_mismatch = [], 0, []
    for name, rows in builder.collect_target(fe_tables).items():
        if name not in ext:
            # legitimate only when the ENDF file curates the whole nuclide
            # away; a present nuclide with a MISSING state name means the
            # state ranking diverged — exactly what this gate must catch
            if _re.sub(r'(g|m\d+)$', '', name) in ext_nuclides:
                state_mismatch.append(name)
            continue
        compared += 1
        e_t = np.array([e for e, _ in rows])
        x_t = np.array([x for _, x in rows])
        e_x = np.array([e for e, _ in ext[name]])
        x_x = np.array([x for _, x in ext[name]])
        x_i = np.interp(e_t, e_x, x_x, left=0.0, right=0.0)
        # a point fails only when it deviates by more than 2% relatively
        # AND by more than 2% of the product's peak absolutely: the s30
        # file legitimately carries up to ~2%-of-peak less than the tables
        # just below the exclusive->MT5 handoff (see the tendl_rp_endf
        # module docstring); a genuine extraction defect (wrong channel,
        # grid smearing, unit slip) exceeds both bounds at once
        rel = np.abs(x_i - x_t) / np.maximum(x_t, 1e-30)
        absd = np.abs(x_i - x_t)
        fail = (x_t > 1e-3) & (rel > 2e-2) & (absd > 2e-2 * x_t.max())
        if fail.any():
            j = int(np.argmax(absd * fail))
            bad.append((name, f'{e_t[j]} MeV: {x_t[j]:.4g} vs {x_i[j]:.4g}'))
    # the gate must not pass vacuously: an empty extraction or a state-set
    # divergence is a failure even though the per-point loop saw nothing
    check(f'rp-from-ENDF-6 extraction reproduces the TENDL-2025 Fe056 tables '
          f'({compared} products compared)',
          not bad and not state_mismatch and compared >= 100,
          str((bad + [('state set', state_mismatch[:5])]) if state_mismatch
              else bad[:3]) if (bad or state_mismatch) else f'{compared} compared')

    mismatch = []
    for curie_tgt, tab, path in endf_targets(src):
        cols, grid, data = assemble_endf(path)
        got = q(con, f'SELECT * FROM {tab} ORDER BY rowid')
        if [tuple(r) for r in got] != [tuple(r) for r in data.tolist()]:
            mismatch.append(tab)
    check('ENDF-6-ingested tables identical to re-extraction', not mismatch,
          str(mismatch[:3]))


def _coverage_diff_rp(con, beam, letter):
    v1p = V1_DIR / f'tendl_{letter}_rp.db'
    if not v1p.exists():
        print(f'  (no v1 tendl_{letter}_rp.db — new library, no diff)')
        return
    v1 = sqlite3.connect(v1p)
    old = {(t, p) for t, p in q(v1, 'SELECT target,product FROM all_reactions')}
    new = {(t, p) for t, p in q(con, 'SELECT target,product FROM all_reactions')}
    old_t, new_t = {t for t, _ in old}, {t for t, _ in new}
    print(f'\n  targets: {len(old_t)} -> {len(new_t)} (+{len(new_t - old_t)}, '
          f'-{len(old_t - new_t)}); dropped: {sorted(old_t - new_t) or "none"}')
    both = old_t & new_t
    lost = [(t, p) for t, p in old if t in both and (t, p) not in new]
    print(f'  products on common targets: '
          f'+{len([x for x in new if x[0] in both and x not in old])} added, '
          f'-{len(lost)} dropped')
    for t, p in lost[:10]:
        print(f'    - {t} -> {p}')


# --- decay ---

DECAY_SCHEMA = {
    'chart': ['isotope', 'name', 'E_level', 'isomer', 'J_pi', 'A', 'Z', 'N',
              'element', 'stable', 'half_life', 'unc_half_life', 'abundance',
              'unc_abundance', 'amu', 'Delta', 'decay_mode'],
    'gammas': ['isotope', 'E_level', 'isomer', 'energy', 'intensity',
               'unc_intensity', 'notes', 'decay_mode', 'unc_norm'],
    'electrons': ['isotope', 'E_level', 'isomer', 'energy', 'intensity',
                  'unc_intensity', 'notes', 'decay_mode'],
    'alphas': ['isotope', 'E_level', 'isomer', 'energy', 'intensity',
               'unc_intensity', 'decay_mode'],
    'beta_minus': ['isotope', 'E_level', 'isomer', 'mean_energy',
                   'intensity', 'unc_intensity', 'endpoint_energy',
                   'decay_mode'],
    'beta_plus': ['isotope', 'E_level', 'isomer', 'mean_energy', 'intensity',
                  'unc_intensity', 'endpoint_energy', 'decay_mode'],
    'NFY': ['parent', 'daughter', 'energy', 'yield', 'unc_yield',
            'cumulative_yield', 'unc_cumulative_yield'],
    'SFY': ['parent', 'daughter', 'yield', 'unc_yield', 'cumulative_yield',
            'unc_cumulative_yield'],
    'coincidences': ['isotope', 'isomer', 'energy', 'coinc_energy',
                     'coinc_type', 'intensity', 'unc_intensity'],
    'elements': ['Z', 'element'],
}


def validate_decay(con):
    for t, cols in DECAY_SCHEMA.items():
        have = [c[1] for c in q(con, f'PRAGMA table_info({t})')]
        check(f'{t} columns', have == cols, str(have))

    names = {r[0] for r in q(con, 'SELECT name FROM chart')}
    dup = q(con, 'SELECT name, COUNT(*) c FROM chart GROUP BY name HAVING c>1')
    check('chart names unique', not dup, str(dup[:3]))

    # decay_mode triples: well-formed, daughters closed, branches <= 1
    bad_triple, missing_daughter, bad_sum = [], [], []
    graph = {}
    for name, dm in q(con, 'SELECT name, decay_mode FROM chart '
                           'WHERE decay_mode IS NOT NULL'):
        tot = 0.0
        for entry in dm.split(','):
            parts = entry.split(':')
            if len(parts) != 3:
                bad_triple.append((name, entry))
                continue
            try:
                br = float(parts[2])
            except ValueError:
                bad_triple.append((name, entry))
                continue
            if not 0.0 < br <= 1.0:
                bad_triple.append((name, entry))
            tot += br
            if parts[1] != 'SFY':
                graph.setdefault(name, []).append(parts[1])
                if parts[1] not in names:
                    missing_daughter.append((name, parts[1]))
        if tot > 1.0 + 1e-5:
            bad_sum.append((name, round(tot, 6)))
    check('decay_mode triples well-formed', not bad_triple, str(bad_triple[:3]))
    check('every decay daughter exists in chart (closure — the DecayChain '
          'guarantee)', not missing_daughter, str(missing_daughter[:5]))
    check('branch fractions sum <= 1 per state', not bad_sum, str(bad_sum[:3]))

    # chains terminate: the daughter graph has no cycles
    WHITE, GRAY, BLACK = 0, 1, 2
    color = dict.fromkeys(graph, WHITE)
    cycle = []

    def dfs(u, stack):
        color[u] = GRAY
        for v in graph.get(u, []):
            if color.get(v, BLACK) == GRAY:
                cycle.append(stack + [v])
                continue
            if color.get(v, BLACK) == WHITE:
                dfs(v, stack + [v])
        color[u] = BLACK
    sys.setrecursionlimit(100000)
    for u in graph:
        if color[u] == WHITE:
            dfs(u, [u])
    check('every decay chain terminates (no cycles)', not cycle,
          str(cycle[:1]))

    # invariants: a NULL-half-life unstable row must not carry radiation
    # that curie's decay-RATE machinery could silently zero — gammas and
    # betas drive activities and count fits; curie instantiates such rows
    # as stable with a warning. Alpha or electron lines alone are allowed
    # (pure spectroscopy lookups: 215Fr's high-spin states are observed
    # through their alpha lines with no measured half-life anywhere), as
    # is a decay_mode string, which the forced-stable read never consults.
    bad_null = []
    for (name,) in q(con, 'SELECT name FROM chart WHERE stable=0 AND '
                          'half_life IS NULL'):
        n_rad = q(
            con,
            'SELECT (SELECT COUNT(*) FROM gammas WHERE isotope||isomer=?) '
            '+ (SELECT COUNT(*) FROM beta_minus WHERE isotope||isomer=?) '
            '+ (SELECT COUNT(*) FROM beta_plus WHERE isotope||isomer=?)',
            (name,) * 3)[0][0]
        if n_rad:
            bad_null.append(name)
    n_null = q(con, 'SELECT COUNT(*) FROM chart WHERE stable=0 AND '
                    'half_life IS NULL')[0][0]
    check(f'unstable NULL-half-life rows ({n_null}) carry no gamma or '
          f'beta rows', not bad_null, str(bad_null[:5]))
    n_amu = q(con, 'SELECT COUNT(*) FROM chart WHERE amu IS NULL')[0][0]
    check('no NULL amu', n_amu == 0, f'{n_amu} rows')

    # cross-table state consistency: every radiation/coincidence state
    # exists in chart with the same E_level
    for t in ('gammas', 'electrons', 'alphas', 'beta_minus', 'beta_plus'):
        bad = q(con, f'SELECT DISTINCT r.isotope, r.isomer FROM {t} r '
                     f'LEFT JOIN chart c ON c.isotope = r.isotope AND '
                     f'c.isomer = r.isomer AND c.E_level = r.E_level '
                     f'WHERE c.name IS NULL')
        check(f'{t} states exist in chart with matching E_level', not bad,
              str(bad[:3]))
    bad = q(con, 'SELECT DISTINCT co.isotope, co.isomer FROM coincidences co '
                 'LEFT JOIN chart c ON c.isotope = co.isotope AND '
                 'c.isomer = co.isomer WHERE c.name IS NULL')
    check('coincidence states exist in chart', not bad, str(bad[:3]))

    # value ranges
    for t in ('gammas', 'electrons', 'alphas', 'beta_minus', 'beta_plus'):
        bad = q(con, f'SELECT COUNT(*) FROM {t} WHERE intensity <= 0 '
                     f'OR intensity > 1000')
        check(f'{t} intensities in (0, 1000] percent', bad[0][0] == 0,
              f'{bad[0][0]} rows')
    grid = [r[0] for r in q(con, 'SELECT DISTINCT energy FROM NFY ORDER BY energy')]
    check('NFY grid byte-exact canonical', grid == [0.0253, 5e5, 2e6, 14e6],
          str(grid))
    bad = q(con, 'SELECT COUNT(*) FROM NFY WHERE ("yield" IS NOT NULL AND '
                 '("yield" < 0 OR "yield" > 1)) OR (cumulative_yield IS NOT '
                 'NULL AND (cumulative_yield < 0 OR cumulative_yield > 2))')
    check('NFY yields in range', bad[0][0] == 0, f'{bad[0][0]} rows')
    fy_names = {r[0] for r in q(con, 'SELECT DISTINCT daughter FROM NFY '
                                     'UNION SELECT daughter FROM SFY '
                                     'UNION SELECT parent FROM NFY '
                                     'UNION SELECT parent FROM SFY')}
    fy_missing = sorted(fy_names - names)
    check('FY parents and daughters exist in chart', not fy_missing,
          str(fy_missing[:5]))

    # spot anchors: certificate lines, half-lives, isomer-routing branches,
    # the preserved-name ledger case, fission yields, coincidence pairs
    def one(sql, args=()):
        r = q(con, sql, args)
        return r[0] if r else None

    n_bp, n_bp_null = q(con, 'SELECT COUNT(*), COUNT(*) - '
                              'COUNT(endpoint_energy) FROM beta_plus')[0]
    check(f'beta_plus endpoint_energy populated ({n_bp} rows)',
          n_bp > 0 and n_bp_null < 0.02 * n_bp, f'{n_bp_null} NULL')

    anchors = [
        ('137CS 661.657 keV = 85.1(2)%',
         "SELECT intensity FROM gammas WHERE isotope='137CS' AND isomer='g' "
         "AND ABS(energy-661.657)<0.05", lambda v: abs(v[0] - 85.1) < 0.3),
        ('60CO 1332.492 keV = 99.9826%',
         "SELECT intensity FROM gammas WHERE isotope='60CO' AND isomer='g' "
         "AND ABS(energy-1332.492)<0.05", lambda v: abs(v[0] - 99.9826) < 0.01),
        ('152EU 344.279 keV = 26.59%',
         "SELECT intensity FROM gammas WHERE isotope='152EU' AND isomer='g' "
         "AND ABS(energy-344.279)<0.05", lambda v: abs(v[0] - 26.59) < 0.2),
        ('152EU 1408.013 keV = 20.87%',
         "SELECT intensity FROM gammas WHERE isotope='152EU' AND isomer='g' "
         "AND ABS(energy-1408.013)<0.05", lambda v: abs(v[0] - 20.87) < 0.2),
        ('60CO half-life 5.2711 y',
         "SELECT half_life FROM chart WHERE name='60COg'",
         lambda v: abs(v[0] - 166344192) < 2e4),
        ('137CS half-life 30.08 y',
         "SELECT half_life FROM chart WHERE name='137CSg'",
         lambda v: abs(v[0] / 3.156e7 - 30.08) < 0.1),
        ('110AGm1 is the 249.83 d state at 117.59 keV (preserved name)',
         "SELECT half_life, E_level FROM chart WHERE name='110AGm1'",
         lambda v: abs(v[0] / 86400 - 249.83) < 0.5
         and abs(v[1] - 0.11759) < 1e-4),
        ('152EUm2 is the 96 min IT state at 147.86 keV (preserved name)',
         "SELECT half_life, E_level FROM chart WHERE name='152EUm2'",
         lambda v: abs(v[0] - 5760) < 120 and abs(v[1] - 0.14786) < 1e-4),
        ('99MO routes 87-88% to 99TCm1',
         "SELECT decay_mode FROM chart WHERE name='99MOg'",
         lambda v: 'Bm1:99TCm1:0.87' in v[0]),
        ('130TE carries its measured two-beta half-life (no 1e-6 s)',
         "SELECT half_life FROM chart WHERE name='130TEg'",
         lambda v: v[0] is None or v[0] > 1e20),
        ('235U thermal NFY 137CS cumulative 6.19%',
         "SELECT cumulative_yield FROM NFY WHERE parent='235Ug' AND "
         "daughter='137CSg' AND energy=0.0253",
         lambda v: abs(v[0] - 0.0619) < 0.002),
        ('22NA beta+ endpoint 545.7 keV',
         "SELECT endpoint_energy FROM beta_plus WHERE isotope='22NA' AND "
         "isomer='g' ORDER BY intensity DESC LIMIT 1",
         lambda v: v[0] is not None and abs(v[0] - 545.7) < 3.0),
        ('60CO coincidence pair 1173&1332 = 99.83(3)',
         "SELECT intensity FROM coincidences WHERE isotope='60CO' AND "
         "isomer='g' AND ABS(energy-1173.228)<0.05 AND "
         "ABS(coinc_energy-1332.492)<0.05 AND coinc_type='g'",
         lambda v: abs(v[0] - 99.834) < 0.1),
        # datasets LiveChart serves without intensities, rebuilt from the
        # ENSDF relative intensities times the normalization record
        ('132LA 464.55 keV = 76(4)% (rebuilt from ENSDF normalization)',
         "SELECT intensity FROM gammas WHERE isotope='132LA' AND "
         "isomer='g' AND ABS(energy-464.55)<0.05",
         lambda v: abs(v[0] - 76.0) < 0.5),
        ('108IN 632.9 keV = 100% (rebuilt from ENSDF normalization)',
         "SELECT intensity FROM gammas WHERE isotope='108IN' AND "
         "isomer='g' AND ABS(energy-632.9)<0.6",
         lambda v: abs(v[0] - 100.0) < 2.0),
    ]
    for label, sql, ok_fn in anchors:
        v = one(sql)
        check(label, v is not None and ok_fn(v), str(v))

    ver = q(con, 'SELECT library, generation FROM _version')
    check('_version stamp', ver == [('decay', 'v2')], str(ver))


def coverage_diff_decay(con):
    v1p = V1_DIR / 'decay.db'
    if not v1p.exists():
        print('  (no v1 file for coverage diff)')
        return
    v1 = sqlite3.connect(v1p)
    old = {r[0] for r in q(v1, 'SELECT name FROM chart')}
    new = {r[0] for r in q(con, 'SELECT name FROM chart')}
    old_iso = {n for n in old if not n.endswith('g')}
    new_iso = {n for n in new if not n.endswith('g')}
    print(f'\n  chart: {len(old)} -> {len(new)} states '
          f'(+{len(new - old)}, -{len(old - new)}); isomers '
          f'{len(old_iso)} -> {len(new_iso)}')
    print(f'  dropped states: {sorted(old - new)[:20]}'
          f'{"..." if len(old - new) > 20 else ""}')

    # half-life moves > 5% on common states (excluding operator-limited)
    moves = []
    v1_hl = dict(q(v1, 'SELECT name, half_life FROM chart'))
    for name, hl in q(con, 'SELECT name, half_life FROM chart'):
        h1 = v1_hl.get(name)
        if h1 and hl and h1 > 0 and abs(hl / h1 - 1) > 0.05:
            moves.append((name, h1, hl))
    print(f'  half-life moves > 5%: {len(moves)}; largest: '
          f'{sorted(moves, key=lambda m: -abs(m[2] / m[1] - 1))[:5]}')

    # branch moves > 0.05 on common states
    def triples(row):
        out = {}
        if row:
            for e in row.split(','):
                p = e.split(':')
                if len(p) == 3:
                    try:
                        out[p[1]] = out.get(p[1], 0.0) + float(p[2])
                    except ValueError:
                        pass
        return out
    v1_dm = dict(q(v1, 'SELECT name, decay_mode FROM chart'))
    br_moves = []
    for name, dm in q(con, 'SELECT name, decay_mode FROM chart'):
        if name not in v1_dm:
            continue
        t1, t2 = triples(v1_dm[name]), triples(dm)
        for d in set(t1) | set(t2):
            if abs(t1.get(d, 0.0) - t2.get(d, 0.0)) > 0.05:
                br_moves.append((name, d, round(t1.get(d, 0.0), 4),
                                 round(t2.get(d, 0.0), 4)))
    print(f'  branch moves > 0.05: {len(br_moves)}; sample: {br_moves[:6]}')

    # gamma-line count changes per common state
    v1_g = dict(q(v1, 'SELECT isotope || isomer, COUNT(*) FROM gammas '
                      'GROUP BY isotope, isomer'))
    v2_g = dict(q(con, 'SELECT isotope || isomer, COUNT(*) FROM gammas '
                       'GROUP BY isotope, isomer'))
    lost_all = [k for k in v1_g if k in new and k not in v2_g]
    print(f'  states with v1 gammas but none in v2: {len(lost_all)} — '
          f'{sorted(lost_all)[:10]}')

    # fission-yield movers (the corrupt-value class)
    v1_y = {(p, d, e): y for p, d, e, y in
            q(v1, 'SELECT parent, daughter, energy, "yield" FROM NFY')}
    big = []
    for p, d, e, y in q(con, 'SELECT parent, daughter, energy, "yield" '
                             'FROM NFY'):
        y1 = v1_y.get((p, d, e))
        if y1 and y and y1 > 1e-4 and (y1 / y > 5 or y / y1 > 5):
            big.append((p, d, e, y1, y))
    print(f'  NFY values moved > 5x (v1 defect class): {len(big)}; sample: '
          f'{[(p, d, e, f"{y1:.3g}", f"{y2:.3g}") for p, d, e, y1, y2 in big[:6]]}')


import functools
VALIDATORS = {
    'decay.db': (validate_decay, coverage_diff_decay),
    'iaea_monitors.db': (validate_iaea_monitors, coverage_diff_iaea),
    'IRDFF.db': (validate_irdff, coverage_diff_irdff),
    'endf.db': (validate_endf, coverage_diff_endf),
    'tendl.db': (validate_tendl, coverage_diff_tendl),
}
for _beam, _l in [('neutron', 'n'), ('proton', 'p'), ('deuteron', 'd'),
                  ('alpha', 'a')]:
    VALIDATORS[f'tendl_{_l}_rp.db'] = (
        functools.partial(_validate_tendl_rp, beam=_beam, letter=_l),
        functools.partial(_coverage_diff_rp, beam=_beam, letter=_l))


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
