"""Build the TENDL residual-production databases (tendl_{n,p,d,a}_rp.db).

Source: the mirrored TENDL-2025 residual tables. Per (target, product):
`rpZZZAAA.tot` is the total production; when isomer-resolved files exist
they are used instead. The `.L##` suffix carries a level INDEX, not an
isomer ordinal: L00 is the ground state and the remaining L-files,
ordered by level index, become m1, m2, ... (an advisory note is printed
when an m-state's header token does not look metastable). A product with
no L-files stores its total under the ground-state name, matching the
previous generation.

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

sys.path.insert(0, str(Path(__file__).parent))
from common import (ARCHIVE_ROOT, ELEMENTS, dedupe_grid, open_fresh,
                    parse_yandf, union_grid_stack, write_version)

DB_NAME = {'neutron': 'tendl_n_rp.db', 'proton': 'tendl_p_rp.db',
           'deuteron': 'tendl_d_rp.db', 'alpha': 'tendl_a_rp.db'}


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
        if lfiles and lfiles[0][0] != 0:
            # Never observed from TALYS (residuals populate the ground
            # state, so an L00 accompanies any L-file) — but if it ever
            # happens, silently dropping the .tot and promoting the lowest
            # isomer to m1 would mislabel states, so stop the build.
            raise SystemExit(f'{tdir}: rp{z:03d}{a:03d} has isomer L-files '
                             f'but no L00 ground file')
        if lfiles:
            chosen, iso_n = [], 0
            for lv, p in lfiles:
                if lv == 0:
                    state = 'g'
                else:
                    iso_n += 1
                    state = f'm{iso_n}'
                chosen.append((state, p))
        else:
            chosen = [('g', p) for _, p in entries]
        for state, p in chosen:
            header, rows = parse_yandf(p.read_text())
            if not rows:
                continue
            name = f'{a}{ELEMENTS[z].upper()}{state}'
            products[name] = dedupe_grid(rows, str(p))
            tok = header.get('nuclide', '')
            if tok and state.startswith('m') and not tok.endswith(('m', 'n')):
                print(f'  note: {p.name} stored as {name} but header '
                      f'token is {tok!r}')
    return products


def assemble(tdir):
    """(curie_target, table_name, cols, grid, data) for one target dir —
    the exact content the build writes, reusable by the validator."""
    m = re.fullmatch(r'([A-Z][a-z]?)(\d{3})(m?)', tdir.name)
    el, a, iso = m.group(1), int(m.group(2)), m.group(3)
    curie_tgt = f'{a}{el.upper()}{iso}'
    products = collect_target(tdir / 'tables' / 'residual')
    if not products:
        return curie_tgt, None, [], None, None
    cols = sorted(products)
    grid, data = union_grid_stack([products[c] for c in cols])
    return curie_tgt, f'{el.upper()}_{a}{iso}', cols, grid, data


def write_target(con, curie_tgt, table, cols, data):
    con.execute(f'CREATE TABLE {table} (energy REAL, ' +
                ', '.join(f'"_{c}" REAL' for c in cols) + ')')
    con.executemany(
        f'INSERT INTO {table} VALUES ({",".join("?" * (len(cols) + 1))})',
        data.tolist())
    con.executemany('INSERT INTO all_reactions VALUES (?,?,?)',
                    [(curie_tgt, c, f'_{c}') for c in cols])
    con.commit()


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
        curie_tgt, table, cols, grid, data = assemble(tdir)
        if table is None:
            print(f'  {curie_tgt}: no residual data, skipped')
            continue
        write_target(con, curie_tgt, table, cols, data)
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
