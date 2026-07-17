"""Build iaea_monitors.db from the mirrored IAEA medical data.

Sources: the 2025 beam-monitor evaluation plus the recommended medical
radionuclide production sets (gamma/positron/therapeutic emitters), all as
per-reaction zips containing tab-separated recommended cross sections
(Pade fit: E [MeV] | sigma [mb] | unc [mb]; yield columns ignored).

Output schema (identical to what curie's Library class reads):
  proton/deuteron/helion/alpha(target, outgoing, product, energy,
                               cross_section, unc_cross_section)
  all_reactions(target, incident, outgoing, product)   -- product is the
      retrieval label and must stay the LAST column
  _sources(target, incident, outgoing, product, family, stem,
           cumulative, mixed_isomer)  -- provenance, not read by curie;
      the two INTEGER flags are 0/1: `cumulative` marks a "(cum)"-labeled
      excitation function, `mixed_isomer` an m+g summed product stored
      under the unsuffixed name
  _version                                              -- provenance stamp

Neutron- and photon-induced entries in the source pages are skipped: the
recommended monitor set covers p/d/3He/alpha, and neutron monitor lookups
resolve through the other libraries.

Rules:
  - reaction identity parsed from the data file's header line, never the
    file name
  - duplicate reactions across families: the monitor evaluation wins; a
    direct evaluation beats its cumulative twin
  - energy grids sorted; exact-duplicate rows collapsed; conflicting
    duplicate abscissae are a build failure unless hand-verified in
    SOURCE_DEFECTS

Usage: python3 iaea_monitors.py
"""
import re
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import ARCHIVE_ROOT, curie_nuclide, curie_target, open_fresh, write_version

FAMILIES = ['monitors', 'gamma', 'positron', 'therapeutic']  # priority order

# Known defects in the IAEA-posted files, each verified by hand. A conflicting
# duplicate abscissa NOT listed here fails the build.
#   rhgp2rh: the last two rows repeat 59.9/60.0 MeV with sigma jumped from
#   76.8 to 99.8 mb - two evaluation segments concatenated at the file tail;
#   the monotone main sequence is kept (keep-first).
SOURCE_DEFECTS = {
    'rhgp2rh': 'keep-first',
}
INCIDENT = {'p': 'p', 'd': 'd', 'a': 'a', '3he': 'h', 'h': 'h',
            'n': 'n', 'g': 'g'}
TABLE = {'p': 'proton', 'd': 'deuteron', 'h': 'helion', 'a': 'alpha'}

_REACTION = re.compile(
    r'^\s*(nat[A-Za-z]{1,2}|\d+[A-Za-z]{1,2})\s*'      # target
    r'\(\s*([^,\s]+)\s*,\s*([^)\s]+)\s*\)\s*'          # (incident, outgoing)
    r'(\d+\s*[mg]?\d*\s*[A-Za-z]{1,2}|\d+[A-Za-z]{1,2}[mg]\d*)\s*$')

# A trailing parenthetical is either a physical-state note ("(gas, 1bar 20C)")
# or the cumulative marker "(cum)"; both are stripped, the latter recorded.
_TRAILING_NOTE = re.compile(r'\s*\(([^()]*)\)\s*$')
# Mixed-isomer products ("87m+gY", "192m1+gIr") are the m+g summed
# excitation function; they map to the UNSUFFIXED product name (matching
# how earlier generations stored the 52Mn/192Ir sums), with the mixed-state
# origin recorded in _sources.
_MIXED_STATE = re.compile(r'(?<=\d)\s*[mg]\d*\+[mg]\d*\s*(?=[A-Za-z])')


def parse_txt(text):
    """(reaction_string, rows[E, sigma, unc]) from a recommended-data txt.

    Energy and cross section are required; the uncertainty column is
    OPTIONAL — several source files leave it blank on the high-energy tail
    (e.g. 125Te(p,n)125I above 25.4 MeV) and the cross section is the
    primary datum. Missing/invalid uncertainties are stored as 0.
    """
    lines = text.splitlines()
    reaction = lines[0].lstrip('﻿').strip()
    rows = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            e, xs = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        try:
            unc = float(parts[2]) if len(parts) > 2 else 0.0
        except ValueError:
            unc = 0.0
        rows.append((e, xs, unc))
    return reaction, rows


def parse_zip(zf):
    """(reaction_string, rows) from a per-reaction zip.

    The .xls member is authoritative: the .txt renders energies with fewer
    decimals, which collides grid points in resonance regions (18O(p,n)18F
    prints 95 duplicate abscissae that the .xls resolves). Falls back to
    the .txt when no .xls is present.
    """
    import xlrd
    xls = [n for n in zf.namelist() if n.lower().endswith('.xls')]
    if not xls:
        txts = [n for n in zf.namelist() if n.lower().endswith('.txt')]
        return parse_txt(zf.read(txts[0]).decode('utf-8', 'replace'))
    sheet = xlrd.open_workbook(file_contents=zf.read(xls[0])).sheet_by_index(0)
    reaction = str(sheet.cell_value(0, 0)).lstrip('﻿').strip()
    rows = []
    for r in range(sheet.nrows):
        e, xs, unc = (sheet.cell_value(r, c) for c in range(3))
        if not (isinstance(e, float) and isinstance(xs, float)):
            continue  # header/blank rows, or no cross section
        if not isinstance(unc, float):
            unc = 0.0  # blank or error cell: sigma is the primary datum
        rows.append((e, xs, unc))
    return reaction, rows


def parse_reaction(s):
    """'natTi(p,x)48V' -> (('natTI', 'p', 'x', '48V'), cumulative, mixed).

    Trailing notes are stripped; '(cum)' sets the cumulative flag; a
    mixed-isomer product ('52m+gMn') maps to the unsuffixed name and sets
    the mixed flag.
    """
    s = s.replace('α', 'a').replace('γ', 'g')
    cumulative = False
    while True:
        m = _TRAILING_NOTE.search(s)
        if not m:
            break
        if 'cum' in m.group(1).lower():
            cumulative = True
        s = s[:m.start()].strip()
    s, n_mixed = _MIXED_STATE.subn('', s)
    mixed = bool(n_mixed)
    m = _REACTION.match(s)
    if not m:
        raise ValueError(f'unparseable reaction string: {s!r}')
    target, inc, outgoing, product = m.groups()
    inc = INCIDENT.get(inc.lower())
    if inc is None:
        raise ValueError(f'unknown incident particle in: {s!r}')
    return ((curie_target(target), inc, outgoing.lower().strip(),
             curie_nuclide(product.replace(' ', ''))), cumulative, mixed)


def collect():
    """Parse every mirrored zip into (reactions, conflicts, skipped).

    `reactions` maps (target, incident, outgoing, product) -> dict(family,
    stem, rows, cumulative, mixed) after family-priority and cumulative-twin
    resolution — the exact content the build writes, reusable by the
    validator for a source-vs-database comparison.
    """
    src = ARCHIVE_ROOT / 'iaea-medical'

    reactions = {}
    skipped, conflicts = [], []

    for family in FAMILIES:
        for zpath in sorted((src / family).glob('*.zip')):
            stem = zpath.stem[:-1]  # drop trailing 't'
            zf = zipfile.ZipFile(zpath)
            if not any(n.lower().endswith(('.xls', '.txt')) for n in zf.namelist()):
                skipped.append((family, stem, 'no data member in zip'))
                continue
            reaction, rows = parse_zip(zf)
            try:
                key, cumulative, mixed = parse_reaction(reaction)
            except ValueError as err:
                skipped.append((family, stem, str(err)))
                continue
            if key[1] in ('n', 'g'):
                skipped.append((family, stem, f'{key[1]}-induced: out of scope'))
                continue
            if not rows:
                skipped.append((family, stem, 'no data rows'))
                continue
            if key in reactions:
                prev = reactions[key]
                # a direct evaluation wins over its cumulative twin (the
                # index key cannot distinguish them); otherwise the earlier
                # (higher-priority) family wins
                if prev['cumulative'] and not cumulative:
                    conflicts.append((key, f'{family} (direct)',
                                      f"{prev['family']} (cumulative twin)"))
                    reactions[key] = {'family': family, 'stem': stem,
                                      'rows': rows, 'cumulative': cumulative,
                                      'mixed': mixed}
                else:
                    why = 'cumulative twin' if cumulative and not prev['cumulative'] \
                        else 'cross-family duplicate'
                    conflicts.append((key, prev['family'], f'{family} ({why})'))
                continue
            reactions[key] = {'family': family, 'stem': stem, 'rows': rows,
                              'cumulative': cumulative, 'mixed': mixed}
    return reactions, conflicts, skipped


def build_grid(key, info, announce=print):
    """Sorted [(E, xs, unc)] for one reaction, applying the exact-duplicate
    collapse and the hand-verified SOURCE_DEFECTS policies; raises on any
    unlisted conflicting duplicate abscissa."""
    target, inc, outgoing, product = key
    grid = {}
    for e, xs, unc in info['rows']:
        if e in grid and grid[e] != (xs, unc):
            if SOURCE_DEFECTS.get(info['stem']) == 'keep-first':
                if announce:
                    announce(f'  known source defect ({info["stem"]}): kept '
                             f'first of conflicting rows at {e} MeV')
                continue
            raise SystemExit(f'conflicting duplicate abscissa {e} MeV in '
                             f'{info["stem"]} ({target}({inc},{outgoing}){product})')
        grid[e] = (xs, unc)
    return [(e, xs, unc) for e, (xs, unc) in sorted(grid.items())]


def main():
    reactions, conflicts, skipped = collect()

    con = open_fresh('iaea_monitors.db')
    for t in TABLE.values():
        con.execute(f'CREATE TABLE {t} (target TEXT, outgoing TEXT, '
                    'product TEXT, energy REAL, cross_section REAL, '
                    'unc_cross_section REAL)')
    con.execute('CREATE TABLE all_reactions (target TEXT, incident TEXT, '
                'outgoing TEXT, product TEXT)')
    con.execute('CREATE TABLE _sources (target TEXT, incident TEXT, '
                'outgoing TEXT, product TEXT, family TEXT, stem TEXT, '
                'cumulative INTEGER, mixed_isomer INTEGER)')

    for key in sorted(reactions):
        target, inc, outgoing, product = key
        info = reactions[key]
        grid = build_grid(key, info)
        con.execute('INSERT INTO all_reactions VALUES (?,?,?,?)',
                    (target, inc, outgoing, product))
        con.execute('INSERT INTO _sources VALUES (?,?,?,?,?,?,?,?)',
                    (target, inc, outgoing, product, info['family'],
                     info['stem'], int(info['cumulative']), int(info['mixed'])))
        con.executemany(
            f'INSERT INTO {TABLE[inc]} VALUES (?,?,?,?,?,?)',
            [(target, outgoing, product, e, xs, unc) for e, xs, unc in grid])

    write_version(con, 'iaea_monitors',
                  'IAEA recommended cross sections (www-nds.iaea.org/medical): '
                  'Monitor Reactions-2025 + gamma/positron/therapeutic emitters',
                  'site as of 2026-07')
    con.commit()

    print(f'{len(reactions)} reactions written')
    by_family = {}
    for info in reactions.values():
        by_family[info['family']] = by_family.get(info['family'], 0) + 1
    print('by family:', by_family)
    mixed = [k for k, i in reactions.items() if i['mixed']]
    cum = [k for k, i in reactions.items() if i['cumulative']]
    print(f'{len(mixed)} m+g summed products stored unsuffixed:',
          ' '.join('%s(%s,%s)%s' % k for k in sorted(mixed)))
    print(f'{len(cum)} cumulative excitation functions:',
          ' '.join('%s(%s,%s)%s' % k for k in sorted(cum)))
    print(f'{len(conflicts)} duplicate/twin resolutions:')
    for key, kept, dropped in conflicts:
        print(f'  {key[0]}({key[1]},{key[2]}){key[3]}: kept {kept}, dropped {dropped}')
    print(f'{len(skipped)} skipped:')
    for s in skipped:
        print('  ', *s)
    con.close()


if __name__ == '__main__':
    main()
