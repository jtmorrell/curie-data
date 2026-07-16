"""Build IRDFF.db from the IAEA IRDFF-II tabulated distribution.

Source: IRDFF-II_TAB.zip -> IRDFF-II.tab.txt, sequential per-reaction
blocks: a header line like 'Tm-169(n,3n)Tm-167' or 'U-235(n,f)' followed
by data lines of ENDF-style floats (E [eV], sigma [barn], absolute
uncertainty [barn], relative uncertainty [%]). The first three columns are
stored; energies and cross sections stay in eV/barn (curie scales on read,
positionally, so column order is load-bearing).

Header conventions: elemental targets ('Ti(n,x)Sc-46') are natural-
composition targets ('natTI'); "(n,n')" is inelastic ('inl'); a trailing
'n' isomer letter is the ENDF second-isomer label ('Ag-110n' = the
117.6 keV, 249.8 d state = m2 in energy-ordered numbering); exclusive
channels may omit the product ('Tm-169(n,3n)'), which is then resolved
against the pinned reaction set.

The reaction set is pinned to the previous generation's 119 (the official
metrology set), with two product names corrected — the old build mangled
the source's isomer letters ('Ag-110n' -> '110NAG', 'Al-26g' -> '26GAL';
now '110AGm2' and '26ALg') — and the build fails if the source no longer
provides any pinned reaction.

Output schema (matching curie's reader):
  _<TARGET>_n_<outgoing>_<PRODUCT>(energy, cross_section, unc_cross_section)
  all_reactions(target, outgoing, product, table_name)
  _version

Usage: python3 irdff.py
"""
import os
import re
import sqlite3
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import ARCHIVE_ROOT, curie_nuclide, open_fresh, write_version

V1_DIR = Path(os.environ.get('CURIE_V1_DIR', Path.home() / '.local/share/curie'))

_HEADER = re.compile(r"^([A-Za-z]{1,2}(?:-\d+)?)\(n,([^)]+)\)"
                     r"\s*([A-Za-z]{1,2}-\d+(?:m\d*|g|n)?)?\s*$")
_ENDF_FLOAT = re.compile(r'(?<=\d)([+-]\d+)$')

# The previous generation corrupted two source isomer letters; the pinned
# reaction list is corrected on read so v2 stores usable names.
V1_NAME_FIXES = {'110NAG': '110AGm2', '26GAL': '26ALg'}


def endf_float(tok):
    """'9.999999-6' -> 9.999999e-6 (ENDF-style exponent without E)."""
    return float(_ENDF_FLOAT.sub(r'E\1', tok))


def nuclide(tok):
    """'Tm-167' / 'In-116m' / 'Ag-110n' -> '167TM' / '116INm' / '110AGm2'."""
    if not tok:
        return ''
    tok = tok.replace('-', '')
    if re.fullmatch(r'[A-Za-z]{1,2}\d+n', tok):
        tok = tok[:-1] + 'm2'   # ENDF second-isomer letter
    return curie_nuclide(tok)


def target_name(tok):
    """'Ti' -> 'natTI' (elemental = natural composition); 'Tm-169' -> '169TM'."""
    return nuclide(tok) if any(c.isdigit() for c in tok) else 'nat' + tok.upper()


def parse_tab(text):
    """{(target, outgoing, product_or_empty): rows[(E, xs, unc)]}.

    A key with an empty product is either a fission entry or a
    product-omitting exclusive-channel header; the caller resolves the
    latter against its pinned reaction list.
    """
    blocks, key = {}, None
    for line in text.splitlines():
        if not line.strip():
            continue
        m = _HEADER.match(line.strip())
        if m:
            tgt, out, prod = m.groups()
            out = out.strip().lower().replace("n'", 'inl')
            key = (target_name(tgt), out, nuclide(prod or ''))
            blocks.setdefault(key, [])
            continue
        if key is None:
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            blocks[key].append((endf_float(parts[0]), endf_float(parts[1]),
                                endf_float(parts[2])))
        except ValueError:
            continue
    return {k: v for k, v in blocks.items() if v}


def main():
    zf = zipfile.ZipFile(ARCHIVE_ROOT / 'irdff-ii' / 'IRDFF-II_TAB.zip')
    blocks = parse_tab(zf.read('IRDFF-II.tab.txt').decode('utf-8', 'replace'))

    v1 = sqlite3.connect(V1_DIR / 'IRDFF.db')
    pinned = sorted(
        (t, o, V1_NAME_FIXES.get(p, p)) for t, o, p in
        v1.execute('SELECT target,outgoing,product FROM all_reactions'))

    # resolve each pinned reaction to a source block: exact key first, then
    # the product-omitting header for the same (target, outgoing) — valid
    # only if it is unambiguous among the pinned reactions
    productless_claims = {}
    resolved, missing = {}, []
    for key in pinned:
        target, outgoing, product = key
        if key in blocks:
            resolved[key] = blocks[key]
        elif product and (target, outgoing, '') in blocks:
            claim = productless_claims.setdefault((target, outgoing), key)
            if claim != key:
                raise SystemExit(f'ambiguous product-omitting source block '
                                 f'{target}(n,{outgoing}): {claim} vs {key}')
            resolved[key] = blocks[(target, outgoing, '')]
        else:
            missing.append(key)
    if missing:
        raise SystemExit(f'source no longer provides pinned reactions: {missing}')

    con = open_fresh('IRDFF.db')
    con.execute('CREATE TABLE all_reactions (target TEXT, outgoing TEXT, '
                'product TEXT, table_name TEXT)')
    for (target, outgoing, product), rows in sorted(resolved.items()):
        table = f'_{target}_n_{outgoing}_{product}'
        con.execute(f'CREATE TABLE {table} (energy REAL, '
                    'cross_section REAL, unc_cross_section REAL)')
        con.execute('INSERT INTO all_reactions VALUES (?,?,?,?)',
                    (target, outgoing, product, table))
        con.executemany(f'INSERT INTO {table} VALUES (?,?,?)', rows)
    write_version(con, 'IRDFF', 'IRDFF-II tabulated (www-nds.iaea.org/IRDFF)',
                  'IRDFF-II final (2020-12-09)')
    con.commit()
    extra = len(blocks) - len(pinned)
    print(f'{len(pinned)} reactions written; {extra} source blocks not in the '
          f'metrology set (elemental/summed entries)')
    con.close()


if __name__ == '__main__':
    main()
