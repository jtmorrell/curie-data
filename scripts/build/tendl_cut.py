"""Derive the curated TENDL target list.

Criteria: targets that are naturally occurring OR have a half-life of at
least one year, with Z <= 100 — applied to ground states and isomeric
targets alike, each isomeric state judged in its own right. An isomeric
target is included when TENDL provides it AND the corresponding chart
state passes the criteria; TALYS isomer-target letters m, n, o track the
curie ordinals m1, m2, m3, so a long-lived second isomer over a short
first one (178HFm2, 192IRm2) is reachable as the 'n' target. 180Ta's
naturally occurring state IS its first isomer, over an 8-hour ground
state. The TENDL-provided enumeration is committed as
tendl_isomer_targets.txt; regenerate it with --probe-isomers, which
checks each candidate state's isomer index page across all four beams.
Half-lives and abundances come from the decay database.

Writes tendl_targets.txt beside this script; run after a decay-database
rebuild or a criteria change to refresh the list.

Usage:
    python3 tendl_cut.py                  # rewrite tendl_targets.txt
    python3 tendl_cut.py --probe-isomers  # re-enumerate TENDL's isomer
                                          # targets first (network)
"""
import re
import sqlite3
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import DECAY_DB

OUT = Path(__file__).parent / 'tendl_targets.txt'
ISOMERS = Path(__file__).parent / 'tendl_isomer_targets.txt'

YEAR = 365.25 * 86400  # seconds


def criteria_cut():
    """-> (admitted ground states, {bare name: admitted isomer ordinals})."""
    dd = sqlite3.connect(DECAY_DB)
    grounds, isomer_ok = set(), {}
    for name, iso, z, stable, hl, ab in dd.execute(
            'SELECT isotope, isomer, Z, stable, half_life, abundance FROM chart'):
        if z > 100:
            continue
        admitted = stable or (ab or 0) > 0 or (hl or 0) >= YEAR
        if not admitted:
            continue
        if iso == 'g':
            grounds.add(name)
        else:
            isomer_ok.setdefault(name, set()).add(int(iso[1:]))
    return grounds, isomer_ok


def _isomer_name(bare, k):
    return bare + ('m' if k == 1 else f'm{k}')


def probe_isomers(grounds, isomer_ok):
    """Enumerate the isomeric targets TENDL provides (any beam's isomer
    index page answers 200) and rewrite the committed enumeration file.
    Candidates are the first-isomer target of every admitted nuclide plus
    the higher-ordinal targets of each isomer that passes the criteria in
    its own right — every name the admission check in main() could
    accept."""
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

    candidates = sorted({n + 'm' for n in grounds | set(isomer_ok)} |
                        {_isomer_name(n, k)
                         for n, ks in isomer_ok.items() for k in ks})
    with ThreadPoolExecutor(12) as pool:
        hits = [t for t, ok in zip(candidates, pool.map(provided, candidates)) if ok]
    ISOMERS.write_text('\n'.join(hits) + '\n')
    print(f'{len(hits)} TENDL-provided isomer targets -> {ISOMERS}')


def main():
    grounds, isomer_ok = criteria_cut()
    if '--probe-isomers' in sys.argv[1:]:
        probe_isomers(grounds, isomer_ok)
    cut = set(grounds)
    if ISOMERS.exists():
        for t in (s.strip() for s in ISOMERS.read_text().splitlines()):
            m = re.fullmatch(r'(.+?)m(\d+)?', t) if t else None
            if m and int(m.group(2) or 1) in isomer_ok.get(m.group(1), ()):
                cut.add(t)
    OUT.write_text('\n'.join(sorted(cut)) + '\n')
    print(f'{len(cut)} targets -> {OUT}')


if __name__ == '__main__':
    main()
