"""Build the TENDL residual-production databases (tendl_{n,p,d,a}_rp.db).

Source: the mirrored TENDL-2025 residual tables. Per (target, product):
`rpZZZAAA.tot` is the total production; when isomer-resolved files exist
(`rpZZZAAA.L00`, `.L01`, ...) they are used instead — L00 is the ground
state and L<i> the i-th isomer, giving products '<A><EL>g', '<A><EL>m1',
... (state names cross-checked against the file header's nuclide token).
A product with no L-files stores its total under the ground-state name,
matching the previous generation.

Output schema per curie's reader: one wide table per target (energy [MeV]
plus one column per product [mb], zero-filled outside each product's
tabulated range on the union grid) and all_reactions(target, product,
label) with the label LAST. The isomer-numbering caveat: L-indices are
TALYS level-ordered isomer numbers; any clash with the decay database's
established names surfaces through the product-resolution check.

Usage: python3 tendl_rp.py <beam>       # neutron|proton|deuteron|alpha
"""
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from common import (ARCHIVE_ROOT, ELEMENTS, open_fresh,
                    parse_yandf, write_version)

DB_NAME = {'neutron': 'tendl_n_rp.db', 'proton': 'tendl_p_rp.db',
           'deuteron': 'tendl_d_rp.db', 'alpha': 'tendl_a_rp.db'}


def target_table(curie_target):
    m = re.fullmatch(r'(\d+)([A-Z]{1,2})(m)?', curie_target)
    a, el, iso = m.groups()
    return f'{el}_{int(a)}{iso or ""}'


def collect_target(tdir):
    """{product_name: rows} for one target's residual directory.

    The .L<n> file suffix is a LEVEL index, not an isomer ordinal
    (rp026052.L42 is 52Fe's first isomer, which happens to sit at level
    42), so state names are assigned by rank: L00 is the ground state and
    the remaining L-files of a product, ordered by level index, are m1,
    m2, ... A product with only a .tot file stores its total under the
    ground-state name.
    """
    files = {}
    for p in sorted(tdir.glob('rp*')):
        m = re.fullmatch(r'rp(\d{3})(\d{3})\.(?:L(\d+)|tot)', p.name)
        if m:
            z, a, level = int(m.group(1)), int(m.group(2)), m.group(3)
            files.setdefault((z, a), []).append(
                (int(level) if level is not None else None, p))
    products = {}
    for (z, a), entries in files.items():
        lfiles = sorted((lv, p) for lv, p in entries if lv is not None)
        if lfiles:
            named = [('g' if lv == 0 else None, lv, p) for lv, p in lfiles]
            iso_n = 0
            chosen = []
            for state, lv, p in named:
                if state is None:
                    iso_n += 1
                    state = f'm{iso_n}'
                chosen.append((state, lv, p))
        else:
            chosen = [('g', None, p) for lv, p in entries]
        for state, level, p in chosen:
            header, rows = parse_yandf(p.read_text())
            if not rows:
                continue
            name = f'{a}{ELEMENTS[z].upper()}{state}'
            products[name] = rows
            tok = header.get('nuclide', '')
            if tok and state.startswith('m') and not tok.endswith(('m', 'n')):
                print(f'  note: {p.name} stored as {name} but header '
                      f'token is {tok!r}')
    return products


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in DB_NAME:
        sys.exit(__doc__.strip())
    beam = sys.argv[1]
    src = ARCHIVE_ROOT / 'tendl-2025' / f'{beam}_files'

    con = open_fresh(DB_NAME[beam])
    con.execute('CREATE TABLE all_reactions (target TEXT, product TEXT, '
                'label TEXT)')
    targets = sorted(t for el in src.iterdir() if el.is_dir()
                     for t in el.iterdir()
                     if (t / 'tables' / 'residual').is_dir())
    n_products = 0
    for tdir in targets:
        nuc = tdir.name  # e.g. Fe056 / Ag110m
        m = re.fullmatch(r'([A-Z][a-z]?)(\d{3})(m?)', nuc)
        curie_tgt = f'{int(m.group(2))}{m.group(1).upper()}{m.group(3)}'
        products = collect_target(tdir / 'tables' / 'residual')
        if not products:
            print(f'  {curie_tgt}: no residual data, skipped')
            continue
        grid = np.unique(np.concatenate([[e for e, _ in r]
                                         for r in products.values()]))
        table = target_table(curie_tgt)
        cols = sorted(products)
        con.execute(f'CREATE TABLE {table} (energy REAL, ' +
                    ', '.join(f'"_{c}" REAL' for c in cols) + ')')
        data = np.column_stack(
            [grid] + [np.interp(grid, [e for e, _ in products[c]],
                                [x for _, x in products[c]], left=0.0,
                                right=0.0) for c in cols])
        con.executemany(
            f'INSERT INTO {table} VALUES ({",".join("?" * (len(cols) + 1))})',
            data.tolist())
        con.executemany('INSERT INTO all_reactions VALUES (?,?,?)',
                        [(curie_tgt, c, f'_{c}') for c in cols])
        con.commit()
        n_products += len(cols)
    write_version(con, DB_NAME[beam].replace('.db', ''),
                  'TENDL-2025 residual-production tables '
                  '(tendl.imperial.ac.uk), curated target set',
                  'TENDL-2025 (TALYS-2.1, 2025-05-14)')
    con.commit()
    print(f'{len(targets)} targets, {n_products} product columns '
          f'-> {DB_NAME[beam]}')
    con.close()


if __name__ == '__main__':
    main()
