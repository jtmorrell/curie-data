"""Add the fission-yield tables to decay.db (NFY + SFY).

Source: the ENDF/B-VIII.1 fission-yield sublibraries — MF8 MT454
(independent yields) and MT459 (cumulative yields) — from the mirrored
per-system files. Rows keep the previous generation's schema with the
additive cumulative columns:

    NFY(parent, daughter, energy [eV], yield, unc_yield,
        cumulative_yield, unc_cumulative_yield)
    SFY(parent, daughter, yield, unc_yield,
        cumulative_yield, unc_cumulative_yield)

Yields are absolute independent (respectively cumulative) fractions per
fission. Energies are stored on the canonical grid ONLY — 0.0253, 5e5,
2e6, 14e6 eV (thermal / fast / 2 MeV / 14 MeV); a system grid point maps
to its class value and the mapping is logged when the source energy
differs (the previous generation's 0.0252 eV stray disappears by
construction). A daughter present only in MT459 carries a NULL
independent yield (flagged); curie's readers select the named columns
they need, so both old and new curie read these tables.

Daughter names resolve against the chart built by decay.py: the FPS field
is the EVALUATION's isomer ordinal, mapped by POSITION onto the chart's
admitted isomers for that nuclide (the chart preserves previous-
generation names, so its ordinals can have gaps — VIII.1's FPS=1 for
73Zn is its 13 ms state, which the chart names 73ZNm2). A daughter
missing from the chart gets a minimal NUBASE2020 row, the same closure
repair the chart itself uses.

Usage: python3 decay_fy.py        # appends to build/decay.db
"""
import re
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import ARCHIVE_ROOT, BUILD_DIR, ELEMENTS
from decay import AMU_KEV, load_nubase

CANONICAL_EV = (0.0253, 5e5, 2e6, 14e6)

fy_report = {'grid_mapped': [], 'cumulative_only': 0, 'independent_only': 0,
             'closure_rows': [], 'fps_position_mapped': [],
             'fps_unresolved': []}


def canonical(e_ev):
    if e_ev < 1.0:
        c = 0.0253
    elif e_ev < 1e6:
        c = 5e5
    elif e_ev < 5e6:
        c = 2e6
    else:
        c = 14e6
    if abs(e_ev - c) > 1e-6 * max(c, 1.0):
        fy_report['grid_mapped'].append((e_ev, c))
    return c


def parent_name(path):
    m = re.fullmatch(r'[ns]fy-(\d+)_([A-Za-z]+)_(\d+)(m\d+)?\.endf',
                     path.name)
    z, el, a, iso = int(m.group(1)), m.group(2), int(m.group(3)), m.group(4)
    assert ELEMENTS[z].lower() == el.lower(), path
    return f'{a}{el.upper()}{iso or "g"}'


class DaughterMap:
    """(ZAFP, FPS) -> chart name, by POSITION among the chart's admitted
    isomers; missing daughters get minimal NUBASE closure rows."""

    def __init__(self, con, nubase):
        self.con = con
        self.nubase = nubase
        self.isomers = {}                 # (A, Z) -> [names by rank order]
        self.known = set()
        for name, a, z, isomer in con.execute(
                'SELECT name, A, Z, isomer FROM chart'):
            self.known.add(name)
            if isomer != 'g':
                self.isomers.setdefault((a, z), []).append(
                    (int(isomer[1:]), name))
        for k in self.isomers:
            self.isomers[k] = [n for _, n in sorted(self.isomers[k])]

    def _closure_row(self, a, z, name, state_idx):
        nb = self.nubase.get((a, z, state_idx))
        nb_gs = self.nubase.get((a, z, 0))
        me = nb_gs['me_kev'] if nb_gs else None
        exc = (nb['exc_kev'] if nb and state_idx else 0.0) or 0.0
        self.con.execute(
            'INSERT INTO chart VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
            (f'{a}{ELEMENTS[z].upper()}', name, exc / 1000.0,
             name[len(f'{a}{ELEMENTS[z].upper()}'):], '', a, z, a - z,
             ELEMENTS[z], 1 if (nb or nb_gs or {}).get('stable') else 0,
             (nb or {}).get('hl_sec'), None, None, None,
             a + me / AMU_KEV if me is not None else None,
             (me / 1000.0 + exc / 1000.0) if me is not None else None, None))
        self.known.add(name)
        fy_report['closure_rows'].append(name)

    def resolve(self, zafp, fps):
        z, a = divmod(int(round(zafp)), 1000)
        base = f'{a}{ELEMENTS[z].upper()}'
        k = int(round(fps))
        if k == 0:
            name = base + 'g'
            if name not in self.known:
                self._closure_row(a, z, name, 0)
            return name
        ranked = self.isomers.get((a, z), [])
        if len(ranked) >= k:
            name = ranked[k - 1]
            if name != f'{base}m{k}':
                fy_report['fps_position_mapped'].append((base, k, name))
            return name
        # no such admitted isomer: NUBASE-backed closure row
        name = f'{base}m{k}'
        if name not in self.known:
            if (a, z, k) in self.nubase:
                self._closure_row(a, z, name, k)
                self.isomers.setdefault((a, z), []).append(name)
            else:
                fy_report['fps_unresolved'].append((base, k))
                self._closure_row(a, z, name, k)   # minimal, hl unknown
                self.isomers.setdefault((a, z), []).append(name)
        return name


def read_yields(path, mt):
    import endf
    sec = endf.Material(str(path)).section_data.get((8, mt))
    return sec['yields'] if sec else []


def merge_system(path, dmap, with_energy):
    """[(parent, daughter, energy|None, y, uy, cy, ucy)] for one file."""
    parent = parent_name(path)
    rows = {}
    for mt, slot in ((454, 0), (459, 1)):
        for entry in read_yields(path, mt):
            e = canonical(entry['E']) if with_energy else None
            for p in entry['products']:
                key = (dmap.resolve(p['ZAFP'], p['FPS']), e)
                y, uy = float(p['Y'][0]), float(p['Y'][1])
                rows.setdefault(key, [None, None, None, None])
                rows[key][2 * slot] = y
                rows[key][2 * slot + 1] = uy
    out = []
    for (daughter, e), (y, uy, cy, ucy) in sorted(rows.items()):
        if y is None:
            fy_report['cumulative_only'] += 1
        if cy is None:
            fy_report['independent_only'] += 1
        out.append((parent, daughter, e, y, uy, cy, ucy))
    return out


def main():
    con = sqlite3.connect(BUILD_DIR / 'decay.db')
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    for t in ('NFY', 'SFY'):
        if t in have:
            con.execute(f'DROP TABLE {t}')
    con.execute('CREATE TABLE NFY (parent TEXT, daughter TEXT, energy REAL, '
                '"yield" REAL, unc_yield REAL, cumulative_yield REAL, '
                'unc_cumulative_yield REAL)')
    con.execute('CREATE TABLE SFY (parent TEXT, daughter TEXT, '
                '"yield" REAL, unc_yield REAL, cumulative_yield REAL, '
                'unc_cumulative_yield REAL)')
    dmap = DaughterMap(con, load_nubase())

    n_rows = s_rows = 0
    nfy_dir = ARCHIVE_ROOT / 'endf-b-viii.1' / 'nfy-version.VIII.1'
    for path in sorted(nfy_dir.glob('nfy-*.endf')):
        rows = merge_system(path, dmap, with_energy=True)
        con.executemany('INSERT INTO NFY VALUES (?,?,?,?,?,?,?)', rows)
        n_rows += len(rows)
    sfy_dir = ARCHIVE_ROOT / 'endf-b-viii.1' / 'sfy-version.VIII.1'
    for path in sorted(sfy_dir.glob('sfy-*.endf')):
        rows = [(p, d, y, uy, cy, ucy)
                for p, d, e, y, uy, cy, ucy in
                merge_system(path, dmap, with_energy=False)]
        con.executemany('INSERT INTO SFY VALUES (?,?,?,?,?,?)', rows)
        s_rows += len(rows)
    con.commit()

    print(f'NFY: {n_rows} rows; SFY: {s_rows} rows')
    print(f'  energies: {[r[0] for r in con.execute("SELECT DISTINCT energy FROM NFY ORDER BY energy")]}')
    if fy_report['grid_mapped']:
        print(f'  grid points mapped to canonical: '
              f'{sorted(set(fy_report["grid_mapped"]))}')
    print(f'  cumulative-only rows (NULL independent yield): '
          f'{fy_report["cumulative_only"]}; independent-only: '
          f'{fy_report["independent_only"]}')
    if fy_report['closure_rows']:
        print(f'  chart closure rows added for FY daughters: '
              f'{len(fy_report["closure_rows"])} — '
              f'{fy_report["closure_rows"][:10]}')
    if fy_report['fps_position_mapped']:
        print(f'  FPS ordinals position-mapped onto chart names: '
              f'{sorted(set(fy_report["fps_position_mapped"]))[:10]}')
    if fy_report['fps_unresolved']:
        print(f'  FPS states with no NUBASE backing (minimal rows): '
              f'{sorted(set(fy_report["fps_unresolved"]))[:10]}')
    con.close()


if __name__ == '__main__':
    main()
