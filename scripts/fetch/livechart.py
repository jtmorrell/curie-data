"""Mirror decay data from the IAEA LiveChart of Nuclides API.

One CSV per (nuclide, endpoint): the levels listing plus decay radiation
per radiation type. The whole-chart ground_states table comes down in a
single request. Already-mirrored files are skipped, so an interrupted
mirror resumes where it left off.

Usage:
    python3 livechart.py --ground-states
    python3 livechart.py NUCLIDE [NUCLIDE ...]     # e.g. 152eu 60co
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from core import fetch

BASE = 'https://nds.iaea.org/relnsd/v1/data?'
RAD_TYPES = ['g', 'bm', 'bp', 'a', 'e', 'x']


def fetch_nuclide(nuc):
    nuc = nuc.lower()
    fetch(f'{BASE}fields=levels&nuclides={nuc}',
          f'livechart/{nuc}', 'levels.csv')
    for rt in RAD_TYPES:
        fetch(f'{BASE}fields=decay_rads&nuclides={nuc}&rad_types={rt}',
              f'livechart/{nuc}', f'decay_rads_{rt}.csv')


if __name__ == '__main__':
    args = sys.argv[1:]
    if '--ground-states' in args:
        fetch(f'{BASE}fields=ground_states&nuclides=all',
              'livechart', 'ground_states.csv')
        args = [a for a in args if a != '--ground-states']
    for nuc in args:
        print(nuc)
        fetch_nuclide(nuc)
