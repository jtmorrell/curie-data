"""Mirror decay data from the IAEA LiveChart of Nuclides API.

One CSV per (nuclide, endpoint): the levels listing plus decay radiation
per radiation type. The whole-chart ground_states table comes down in a
single request; `--all` then iterates every nuclide it lists. A nuclide
with no data for an endpoint returns a tiny numeric error-code body, which
is archived as-is (it is the documented "nothing here" marker). Verified
already-mirrored files are skipped, so an interrupted mirror resumes.

Usage:
    python3 livechart.py --ground-states
    python3 livechart.py NUCLIDE [NUCLIDE ...]     # e.g. 152eu 60co
    python3 livechart.py --all                     # full mirror (hours)
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from core import ARCHIVE_ROOT, fetch

BASE = 'https://nds.iaea.org/relnsd/v1/data?'
RAD_TYPES = ['g', 'bm', 'bp', 'a', 'e', 'x']


def all_nuclides():
    """Nuclide query tokens (e.g. '135xe') from the ground_states table."""
    gs = fetch(f'{BASE}fields=ground_states&nuclides=all',
               'livechart', 'ground_states.csv')
    tokens = []
    with open(gs, newline='') as f:
        for row in csv.DictReader(f):
            try:
                a = int(row['z']) + int(row['n'])
            except (KeyError, ValueError):
                continue
            tokens.append(f"{a}{row['symbol'].strip().lower()}")
    return tokens


def fetch_nuclide(nuc):
    nuc = nuc.lower()
    fetch(f'{BASE}fields=levels&nuclides={nuc}',
          f'livechart/{nuc}', 'levels.csv')
    for rt in RAD_TYPES:
        fetch(f'{BASE}fields=decay_rads&nuclides={nuc}&rad_types={rt}',
              f'livechart/{nuc}', f'decay_rads_{rt}.csv')


if __name__ == '__main__':
    args = sys.argv[1:]
    known_flags = {'--ground-states', '--all'}
    if not args or '-h' in args or '--help' in args or \
            any(a.startswith('-') for a in args if a not in known_flags):
        sys.exit(__doc__.strip())
    if '--all' in args:
        nucs = all_nuclides()
        print(f'{len(nucs)} nuclides to mirror')
        failures = []
        for i, nuc in enumerate(nucs):
            try:
                fetch_nuclide(nuc)
            except Exception as e:
                failures.append((nuc, repr(e)))
                print(f'{nuc}: FAILED {e}', flush=True)
            if (i + 1) % 100 == 0:
                print(f'... {i + 1}/{len(nucs)}', flush=True)
        print(f'done: {len(nucs) - len(failures)} mirrored, {len(failures)} failures')
        for f in failures:
            print('  FAILED:', *f)
        sys.exit(0)
    if '--ground-states' in args:
        fetch(f'{BASE}fields=ground_states&nuclides=all',
              'livechart', 'ground_states.csv')
        args = [a for a in args if a != '--ground-states']
    for nuc in args:
        print(nuc)
        fetch_nuclide(nuc)
