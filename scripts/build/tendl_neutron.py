"""Build tendl.db (neutron exclusive channels) from the TENDL-2025 tables.

Source: the mirrored `tables/xs/` files per target. `xsNPDTHA.tot` encodes
the emitted-particle multiplicities in its six digits (n, p, d, t, h, a) —
`xs200000` is (n,2n), `xs011002` is (n,pd2a) — which is also how the
column label and outgoing string are formed ('2n', 'pd2a'); the special
files map to the summary and production columns (total -> _tot_, elastic
-> _el_, nonelastic -> _non_, fission -> _f_, xs000000 -> _g_ capture,
<x>prod -> _x_<x>). `all.tot` and `reaction.tot` are redundant sums and
are skipped. Residual products for emission channels are computed from
the particle balance, stored as ground-state names.

Channels are stored on the union of their grids (the production files
use threshold-truncated subsets of TALYS's master grid), zero-filled
outside each channel's range; energy [MeV], cross sections [mb], the
units curie reads without scaling. Value-identical duplicate rows (a
round-off defect the source still carries) are collapsed; conflicting
duplicate abscissae are a build failure (curie's reader would collapse
them silently).

Usage: python3 tendl_neutron.py
"""
import re
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (ARCHIVE_ROOT, ELEMENTS, open_fresh, parse_yandf,
                    write_version)

PARTICLE_AZ = {'n': (1, 0), 'p': (1, 1), 'd': (2, 1), 't': (3, 1),
               'h': (3, 2), 'a': (4, 2)}
SPECIALS = {'total.tot': ('_tot_', 'tot'), 'elastic.tot': ('_el_', 'el'),
            'nonelastic.tot': ('_non_', 'non'), 'fission.tot': ('_f_', 'f'),
            'nprod.tot': ('_x_n', 'x'), 'gprod.tot': ('_x_g', 'x'),
            'pprod.tot': ('_x_p', 'x'), 'dprod.tot': ('_x_d', 'x'),
            'tprod.tot': ('_x_t', 'x'), 'hprod.tot': ('_x_h', 'x'),
            'aprod.tot': ('_x_a', 'x')}
SKIP = {'all.tot', 'reaction.tot'}


def channel_string(digits):
    """'011002' -> 'pd2a'; '000000' -> '' (capture)."""
    return ''.join((str(m) if m > 1 else '') + letter
                   for m, letter in zip(map(int, digits), 'npdtha') if m > 0)


def residual(a, z, digits):
    pa, pz = a + 1, z
    for m, letter in zip(map(int, digits), 'npdtha'):
        da, dz = PARTICLE_AZ[letter]
        pa, pz = pa - m * da, pz - m * dz
    return f'{pa}{ELEMENTS[pz].upper()}' if 0 < pz < len(ELEMENTS) else ''


def main():
    src = ARCHIVE_ROOT / 'tendl-2025' / 'neutron_files'
    con = open_fresh('tendl.db')
    con.execute('CREATE TABLE all_reactions (target TEXT, outgoing TEXT, '
                'product TEXT, label TEXT)')
    targets = sorted(t for el in src.iterdir() if el.is_dir()
                     for t in el.iterdir() if (t / 'tables' / 'xs').is_dir())
    unknown = set()
    for tdir in targets:
        m = re.fullmatch(r'([A-Z][a-z]?)(\d{3})(m?)', tdir.name)
        el, a, iso = m.group(1), int(m.group(2)), m.group(3)
        curie_tgt = f'{a}{el.upper()}{iso}'
        z = ELEMENTS.index(el)
        cols, grid = [], None
        for p in sorted((tdir / 'tables' / 'xs').iterdir()):
            if p.name in SKIP:
                continue
            xm = re.fullmatch(r'xs(\d{6})\.tot', p.name)
            if p.name in SPECIALS:
                label, out = SPECIALS[p.name]
                prod = label[-1] if out == 'x' else ''
            elif xm:
                digits = xm.group(1)
                out = channel_string(digits) or 'g'
                label = f'_{out}_'
                prod = residual(a, z, digits)
            else:
                unknown.add(p.name)
                continue
            _, rows = parse_yandf(p.read_text())
            if not rows:
                continue
            seen = {}
            for e, x in rows:
                if e in seen and seen[e] != x:
                    raise SystemExit(f'conflicting duplicate abscissa {e} in '
                                     f'{tdir.name}/{p.name}')
                seen[e] = x
            cols.append((label, out, prod, sorted(seen.items())))
        if not cols:
            print(f'  {curie_tgt}: no xs data, skipped')
            continue
        grid = np.unique(np.concatenate([[e for e, _ in r]
                                         for *_, r in cols]))
        table = f'{el.upper()}_{a}{iso}'
        con.execute(f'CREATE TABLE {table} (energy REAL, ' +
                    ', '.join(f'"{lab}" REAL' for lab, *_ in cols) + ')')
        data = np.column_stack(
            [grid] + [np.interp(grid, [e for e, _ in r], [x for _, x in r],
                                left=0.0, right=0.0) for *_, r in cols])
        con.executemany(
            f'INSERT INTO {table} VALUES ({",".join("?" * (len(cols) + 1))})',
            data.tolist())
        con.executemany('INSERT INTO all_reactions VALUES (?,?,?,?)',
                        [(curie_tgt, out, prod, lab) for lab, out, prod, _ in cols])
        con.commit()
    write_version(con, 'tendl',
                  'TENDL-2025 exclusive-channel tables '
                  '(tendl.imperial.ac.uk), curated target set',
                  'TENDL-2025 (TALYS-2.1, 2025-05-14)')
    con.commit()
    print(f'{len(targets)} targets -> tendl.db')
    if unknown:
        print('unrecognized xs files (not stored):', sorted(unknown))
    con.close()


if __name__ == '__main__':
    main()
