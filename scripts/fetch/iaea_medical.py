"""Mirror the IAEA recommended charged-particle cross-section data.

Two families, one site (www-nds.iaea.org/medical/):
  - beam-monitor reactions (monitor_reactions.html, the 2025 evaluation), and
  - medical radionuclide production data (gamma_emitters.html,
    positron_emitters.html, therapeutic.html).

Each reaction has a frameset page `<stem>0.html` and a numerical-download
archive `<stem>t.zip` containing the recommended (Pade-fit) cross section
as tab-separated text. The zip's reaction identity comes from the first
line of the .txt inside it, not from the file name.

Usage: python3 iaea_medical.py
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from core import fetch, fetch_text

BASE = 'https://www-nds.iaea.org/medical/'
INDEX_PAGES = {
    'monitors': 'monitor_reactions.html',
    'gamma': 'gamma_emitters.html',
    'positron': 'positron_emitters.html',
    'therapeutic': 'therapeutic.html',
}


def reaction_stems(index_html):
    """Per-reaction page stems referenced by an index page (`<stem>0.html`).

    The stem pattern is strictly alphanumeric, which also excludes the
    archived-2007 monitor list (it lives under a subdirectory path).
    """
    return sorted(set(re.findall(r'href="([a-z0-9]+)0\.html"', index_html)))


def main():
    failures = []
    seen = {}
    for family, page in INDEX_PAGES.items():
        html = fetch_text(BASE + page, 'iaea-medical', page)
        stems = reaction_stems(html)
        print(f'{family}: {len(stems)} reaction pages')
        for stem in stems:
            if stem in seen:
                print(f'  {stem}: already mirrored (also in {seen[stem]})')
                continue
            try:
                fetch(f'{BASE}{stem}t.zip', f'iaea-medical/{family}')
            except Exception as e:
                # not added to `seen`: a stem that fails here is retried if
                # another family's index also lists it
                failures.append((family, stem, repr(e)))
                print(f'  {stem}t.zip: FAILED {e}')
                continue
            seen[stem] = family
    print(f'\n{len(seen)} unique reactions mirrored, {len(failures)} failures')
    for f in failures:
        print('  MISSING:', *f)


if __name__ == '__main__':
    main()
