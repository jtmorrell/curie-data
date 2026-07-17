"""Derive the curated TENDL target list.

Criteria: ground states that are naturally occurring OR have a half-life
of at least one day, with Z <= 100, unioned with every target the previous
generation shipped (including its isomeric targets) and with the isomeric
targets TENDL itself provides wherever the criteria admit the ground state
or the isomer in its own right (180Ta's naturally occurring state IS the
isomer, over an 8-hour ground state). The TENDL-provided enumeration is
committed as tendl_isomer_targets.txt; regenerate it with --probe-isomers,
which checks each admitted ground state's isomer index page across all
four beams. Half-lives and abundances come from the decay database.

Writes tendl_targets.txt beside this script; run after a decay-database
rebuild to refresh the list.

Usage:
    python3 tendl_cut.py                  # rewrite tendl_targets.txt
    python3 tendl_cut.py --probe-isomers  # re-enumerate TENDL's isomer
                                          # targets first (network)
"""
import os
import sqlite3
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

V1_DIR = Path(os.environ.get('CURIE_V1_DIR', Path.home() / '.local/share/curie'))
OUT = Path(__file__).parent / 'tendl_targets.txt'
ISOMERS = Path(__file__).parent / 'tendl_isomer_targets.txt'


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


def criteria_cut():
    """-> (admitted ground states, bare names whose isomer itself passes)."""
    dd = sqlite3.connect(V1_DIR / 'decay.db')
    grounds, isomer_ok = set(), set()
    for name, iso, z, stable, hl, ab in dd.execute(
            'SELECT isotope, isomer, Z, stable, half_life, abundance FROM chart'):
        if z > 100:
            continue
        admitted = stable or (ab or 0) > 0 or (hl or 0) >= 86400
        if not admitted:
            continue
        (grounds if iso == 'g' else isomer_ok).add(name)
    return grounds, isomer_ok


def probe_isomers(bare_names):
    """Enumerate the isomeric targets TENDL provides for these bare
    nuclide names (any beam's isomer index page answers 200) and rewrite
    the committed enumeration file. Candidates must span every name the
    admission check in main() would accept — the grounds AND the nuclides
    admitted through their isomer (180Ta), or regeneration would silently
    drop the latter."""
    sys.path.insert(0, str(Path(__file__).parents[1] / 'fetch'))
    from tendl import BEAMS, http_get, index_pages

    def provided(iso_target):
        for beam in BEAMS:
            try:
                http_get(index_pages(beam, iso_target)[0])
                return True
            except FileNotFoundError:
                continue        # 404: this beam genuinely lacks the target
        return False            # a transient network error propagates instead

    candidates = sorted(n + 'm' for n in bare_names)
    with ThreadPoolExecutor(12) as pool:
        hits = [t for t, ok in zip(candidates, pool.map(provided, candidates)) if ok]
    ISOMERS.write_text('\n'.join(hits) + '\n')
    print(f'{len(hits)} TENDL-provided isomer targets -> {ISOMERS}')


def main():
    grounds, isomer_ok = criteria_cut()
    if '--probe-isomers' in sys.argv[1:]:
        probe_isomers(grounds | isomer_ok)
    cut = grounds | v1_targets()
    if ISOMERS.exists():
        cut |= {t for t in (s.strip() for s in ISOMERS.read_text().splitlines())
                if t and (t.rstrip('m') in grounds or t.rstrip('m') in isomer_ok)}
    OUT.write_text('\n'.join(sorted(cut)) + '\n')
    print(f'{len(cut)} targets -> {OUT}')


if __name__ == '__main__':
    main()
