"""Derive the curated TENDL target list.

Criteria: ground states that are naturally occurring OR have a half-life
of at least one day, with Z <= 100, unioned with every target the previous
generation shipped (including its isomeric targets). Half-lives and
abundances come from the decay database.

Writes tendl_targets.txt beside this script; run after a decay-database
rebuild to refresh the list.

Usage: python3 tendl_cut.py
"""
import os
import sqlite3
from pathlib import Path

V1_DIR = Path(os.environ.get('CURIE_V1_DIR', Path.home() / '.local/share/curie'))
OUT = Path(__file__).parent / 'tendl_targets.txt'


def v1_targets():
    out = set()
    for db in ['tendl', 'tendl_n_rp', 'tendl_p_rp', 'tendl_d_rp']:
        con = sqlite3.connect(V1_DIR / f'{db}.db')
        for (t,) in con.execute("SELECT name FROM sqlite_master WHERE type='table'"):
            if t == 'all_reactions':
                continue
            el, a = t.split('_')
            out.add(a.rstrip('m') + el + ('m' if a.endswith('m') else ''))
    return out


def main():
    dd = sqlite3.connect(V1_DIR / 'decay.db')
    cut = set()
    for name, iso, z, stable, hl, ab in dd.execute(
            'SELECT isotope, isomer, Z, stable, half_life, abundance FROM chart'):
        if z > 100 or iso != 'g':
            continue
        if stable or (ab or 0) > 0 or (hl or 0) >= 86400:
            cut.add(name)
    cut |= v1_targets()
    OUT.write_text('\n'.join(sorted(cut)) + '\n')
    print(f'{len(cut)} targets -> {OUT}')


if __name__ == '__main__':
    main()
