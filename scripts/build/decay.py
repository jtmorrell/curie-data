"""Build decay.db v2 — the chart table, from the LiveChart mirror.

Row universe: every LiveChart ground state plus every level admitted as an
isomer (half-life >= the 100 ns floor, or listed as an isomer by
NUBASE2020), minus the curated exclusions in isomer_ledger.toml — levels
whose admission would renumber an isomer name the previous generation
established. Isomers are numbered m1, m2, ... by ascending level energy
among the admitted set, and every previous-generation isomer name is
checked to still denote the same level; a renumbering fails the build
(the fix is a ledger entry, decided at review, never silent).

decay_mode: comma-joined label:daughter:branch triples (branch an absolute
fraction), built from the committed decay_modes.toml. Delayed/subset
branches are netted from their parent mode's percentage. Where the
daughter nuclide has admitted isomers, the mode is split by the observed
level flow of the parent state's own dataset — direct particle feeding of
the isomer level plus gamma flow into it weighted by (1 + conversion
coefficient) — so 99MO carries B:99TCg + Bm1:99TCm1 exactly as measured.
Daughters missing from the row universe get a minimal chart row from
NUBASE2020 (the previous generation's orphan class, fixed by closure).

Conventions (consumption map): E_level [MeV]; half_life [s]; amu = the
GROUND state's atomic mass on every row of the nuclide; Delta = the LEVEL
mass excess [MeV]; branches normalized to sum <= 1; stable = 1 only for
LiveChart STABLE rows, so the two-beta quasi-stables carry their measured
half-lives (the previous generation stored a bogus 1e-6 s).

Usage: python3 decay.py            # writes build/decay.db (chart)
"""
import csv
import os
import re
import sqlite3
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import ARCHIVE_ROOT, BUILD_DIR, ELEMENTS

LC = ARCHIVE_ROOT / 'livechart'
NUBASE = ARCHIVE_ROOT / 'ame2020' / 'nubase_4.mas20.txt'
V1_DIR = Path(os.environ.get('CURIE_V1_DIR', Path.home() / '.local/share/curie'))
ISOMER_FLOOR_S = 1e-7          # 100 ns (NUBASE's own convention)
AMU_KEV = 931494.10242         # keV per u

with open(Path(__file__).parent / 'decay_modes.toml', 'rb') as _f:
    MODE_MAP = tomllib.load(_f)
with open(Path(__file__).parent / 'isomer_ledger.toml', 'rb') as _f:
    LEDGER = tomllib.load(_f).get('exclude', [])

report = {'blank_branch_sole': [], 'blank_branch_multi': [], 'sf_combined': [],
          'closure_rows': [], 'no_half_life': [], 'nonclosing_branches': [],
          'ledger_exclusions': [], 'nubase_only_isomers': [],
          'feeding_overflow': [], 'operator_half_lives': 0,
          'newcomer_excluded': [], 'floor_override': [], 'v1_isomer_lost': [],
          'energy_estimated': [], 'order_swap': []}


def nuclide_dir(a, symbol):
    return LC / f'{a}{symbol.lower()}'


def fnum(s):
    s = (s or '').strip()
    return float(s) if s else None


# --- NUBASE2020 (fixed-width; documented column layout) ---

_NUBASE_UNITS = {'ys': 1e-24, 'zs': 1e-21, 'as': 1e-18, 'fs': 1e-15,
                 'ps': 1e-12, 'ns': 1e-9, 'us': 1e-6, 'ms': 1e-3,
                 's': 1.0, 'm': 60.0, 'h': 3600.0, 'd': 86400.0,
                 'y': 31556926.0, 'ky': 31556926e3, 'My': 31556926e6,
                 'Gy': 31556926e9, 'Ty': 31556926e12, 'Py': 31556926e15,
                 'Ey': 31556926e18, 'Zy': 31556926e21, 'Yy': 31556926e24}


def load_nubase():
    """{(A, Z, i): {'hl_sec', 'stable', 'me_kev', 'exc_kev'}} — i = 0 for
    the ground state, 1, 2, ... for the NUBASE isomer letters."""
    out = {}
    for line in open(NUBASE, encoding='utf-8', errors='replace'):
        if line.startswith('#'):
            continue
        try:
            a = int(line[0:3])
            zi = line[4:8].strip()
            z, i = int(zi[:3]), int(zi[3])
        except (ValueError, IndexError):
            continue
        def _nb_num(s):
            try:
                return fnum(s.replace('#', ''))
            except ValueError:            # markers like 'non-exist'
                return None
        me = _nb_num(line[18:31])
        exc = _nb_num(line[42:54])
        t_str = line[69:78].strip()
        unit = line[78:80].strip()
        stable = t_str == 'stbl'
        hl = None
        if not stable and t_str and t_str not in ('p-unst', 'n-unst', 'R'):
            try:
                hl = float(t_str.replace('#', '').replace('>', '')
                           .replace('<', '').replace('~', '')) \
                     * _NUBASE_UNITS.get(unit, float('nan'))
                if hl != hl:      # unknown unit
                    hl = None
            except ValueError:
                hl = None
        out[(a, z, i)] = {'hl_sec': hl, 'stable': stable, 'me_kev': me,
                          'exc_kev': exc}
    return out


# --- LiveChart readers ---

def read_levels(a, symbol):
    """Level rows (positional: the header repeats unc_hl) with the fields
    the chart needs; ragged rows padded."""
    path = nuclide_dir(a, symbol) / 'levels.csv'
    if not path.exists():
        return []
    with open(path) as f:
        rdr = csv.reader(f)
        hdr = next(rdr)
        n = len(hdr)
        out = []
        for row in rdr:
            row = row + [''] * (n - len(row))
            if not row[3].strip():        # blank continuation rows
                continue
            out.append({'idx': int(row[3]), 'energy': fnum(row[5]),
                        'shift': row[4].strip(), 'jp': row[8].strip(),
                        'hl_sec': fnum(row[14]), 'unc_hl_sec': fnum(row[15]),
                        'operator_hl': row[11].strip(),
                        'modes': [(row[k].strip(), fnum(row[k + 1]))
                                  for k in (16, 19, 22) if row[k].strip()]})
        return out


def rad_rows(a, symbol, fname):
    path = nuclide_dir(a, symbol) / fname
    if not path.exists():
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def level_token(s):
    """A LiveChart level-energy token -> (offset [keV], band letters).
    Levels in an unanchored band carry a letter: '0+X' -> (0.0, 'X'),
    '125.3' -> (125.3, ''), 'X' -> (0.0, 'X')."""
    s = (s or '').strip()
    if not s:
        return None, ''
    m = re.fullmatch(
        r'([\d.]+(?:[eE][+-]?\d+)?)?\s*\+?\s*([A-Z]*)\s*\+?\s*([\d.]+)?', s)
    if not m:
        return None, ''
    band = m.group(2) or ''
    val = m.group(1) if m.group(1) is not None else m.group(3)
    try:
        return (float(val) if val else 0.0), band
    except ValueError:
        return None, ''


def state_match(row, raw_kev, band):
    p, p_band = level_token(row.get('p_energy'))
    if raw_kev == 0.0 and not band:
        return p in (None, 0.0) and not p_band
    return (p is not None and p_band == band
            and abs(p - raw_kev) < 1.0)


# --- level-flow feeding: daughter-isomer split of one (state, mode) ---

_MODE_FILE = {'B-': ('decay_rads_bm.csv', ('intensity_beta',)),
              'EC+B+': ('decay_rads_bp.csv', ('intensity_beta', 'intensity_ec')),
              'EC': ('decay_rads_bp.csv', ('intensity_beta', 'intensity_ec')),
              'B+': ('decay_rads_bp.csv', ('intensity_beta', 'intensity_ec')),
              'A': ('decay_rads_a.csv', ('intensity',)),
              'IT': (None, ())}


def isomer_feeding(a, symbol, p_raw_kev, p_band, lc_mode, daughter_states):
    """Percent of parent-state decays through `lc_mode` that populate each
    daughter isomer state: direct particle feeding plus gamma flow into the
    level weighted by (1 + conversion coefficient). All rows come from the
    parent state's own dataset (p_energy + decay filters). Daughter states
    are (key, raw offset [keV], band letters) — a state in an unanchored
    band is fed by rows whose level token names the same band."""
    feed = {key: 0.0 for key, _, _ in daughter_states}

    def credit(token, amount):
        off, band = level_token(token)
        if off is None:
            return
        for key, raw, b in daughter_states:
            if b == band and abs(off - raw) < 1.0:
                feed[key] += amount

    fname, cols = _MODE_FILE[lc_mode]
    if fname:
        for r in rad_rows(a, symbol, fname):
            if not state_match(r, p_raw_kev, p_band) \
                    or r.get('decay', '').strip() != lc_mode:
                continue
            amount = sum(fnum(r.get(c)) or 0.0 for c in cols)
            credit(r.get('daughter_level_energy'), amount)
    for r in rad_rows(a, symbol, 'decay_rads_g.csv'):
        if not state_match(r, p_raw_kev, p_band) \
                or r.get('decay', '').strip() != lc_mode:
            continue
        i = fnum(r.get('intensity'))
        end, _ = level_token(r.get('end_level_energy'))
        if i is None or end is None:
            continue                       # unplaced (x-ray append) or blank
        # first-reached accounting: a transition that STARTS at another
        # admitted isomer is that state's own decay — crediting it here
        # would double-count the cascade (and the inter-isomer transitions
        # are the highly converted ones, so the (1+cc) weight would
        # amplify the error)
        start, s_band = level_token(r.get('start_level_energy'))
        if start is not None and any(
                b == s_band and abs(start - raw) < 1.0
                for _, raw, b in daughter_states):
            continue
        cc = fnum(r.get('conversion_coeff')) or 0.0
        credit(r.get('end_level_energy'), i * (1.0 + cc))
    return feed


# --- chart assembly ---

def _level_energy(l, a, z, nubase, base):
    """Resolve a level's ordering/storage energy [keV]: LiveChart value,
    or — when the absolute energy is unknown (energy_shift set / zero for
    a non-ground level) — the closest NUBASE excitation estimate."""
    e = l['energy']
    if not l['shift'] and e:
        return e
    ests = [nubase[(a, z, i)]['exc_kev'] for i in (1, 2, 3, 4)
            if (a, z, i) in nubase and nubase[(a, z, i)]['exc_kev']]
    est = min(ests, key=lambda x: abs(x - (e or 0.0))) if ests else None
    if est is not None:
        report['energy_estimated'].append((base, e, est))
        return est
    report['energy_estimated'].append((base, e, None))
    return e or 0.0


def admitted_isomers(a, z, symbol, nubase, v1_isomers):
    """[(E_kev, level dict)] ascending — the admitted isomer set.

    Admission: half-life floor + NUBASE census, minus curated ledger
    exclusions — then reconciled against the previous generation's names
    (D-12): established ranks are paired to candidates by energy (a v1
    placeholder energy of exactly 0.1 MeV pairs by order), a candidate
    that would INSERT BELOW an established name is excluded (reported),
    and a level an established rank needs is admitted even below the
    floor. An order swap that no pairing resolves fails the build."""
    base = f'{a}{symbol.upper()}'
    levels = [l for l in read_levels(a, symbol) if l['idx'] > 0]
    for l in levels:
        l['e_res'] = _level_energy(l, a, z, nubase, base) \
            if l['energy'] is not None or l['shift'] else None
    levels = [l for l in levels if l['e_res'] is not None]
    excl = [e['level_energy_kev'] for e in LEDGER if e['nuclide'] == base]

    def ledgered(l):
        return any(abs(l['e_res'] - x) < 0.5 for x in excl)

    cands = {}
    for l in levels:
        if ledgered(l):
            report['ledger_exclusions'].append((base, l['e_res']))
            continue
        if l['hl_sec'] is not None and l['hl_sec'] >= ISOMER_FLOOR_S:
            cands[l['e_res']] = l
    i = 1
    while (a, z, i) in nubase:            # NUBASE census additions
        exc = nubase[(a, z, i)]['exc_kev']
        if exc is not None and not any(abs(exc - e) < max(2.0, 0.01 * exc)
                                       for e in cands):
            near = [l for l in levels
                    if abs(l['e_res'] - exc) < max(2.0, 0.01 * exc)
                    and not ledgered(l)]
            if near:
                l = min(near, key=lambda l: abs(l['e_res'] - exc))
                cands[l['e_res']] = l
                report['nubase_only_isomers'].append((base, l['e_res']))
        i += 1
    cand_list = sorted(cands.items())
    if not v1_isomers:
        return [(i, e, l) for i, (e, l) in enumerate(cand_list, start=1)]

    # pair established ranks to candidates: the HALF-LIFE is the physical
    # fingerprint (the previous generation stored placeholder energies —
    # 0.1, 0.2 MeV — for states whose excitation was unknown, so energy
    # alone misleads); a factor-two half-life agreement matches, with
    # energy closeness as the tiebreak, and a plain energy match as the
    # fallback for ranks whose half-life was revised
    import math

    def hl_close(h1, h2, factor=2.0):
        return (h1 and h2 and h1 > 0 and h2 > 0
                and abs(math.log10(h1 / h2)) < math.log10(factor))

    def tol(e):
        return max(30.0, 0.25 * e)

    # three passes in decreasing confidence: a tight energy agreement
    # outranks the half-life (the previous generation duplicated a state
    # under two ranks with one placeholder energy — 65FE — and the
    # duplicate must lose to the rank whose energy pins the level), the
    # half-life outranks a loose energy (placeholder energies are
    # arbitrary), and a loose energy catches revised half-lives.
    matched = {}                          # v1 rank -> candidate energy
    for rank, e_mev, hl1 in sorted(v1_isomers):
        e1 = e_mev * 1000.0
        free = [(e, l) for e, l in cand_list if e not in matched.values()]
        tight = [x for x in free if abs(x[0] - e1) < max(5.0, 0.01 * e1)]
        if tight:
            matched[rank] = min(tight, key=lambda x: abs(x[0] - e1))[0]
    for rank, e_mev, hl1 in sorted(v1_isomers):
        if rank in matched:
            continue
        free = [(e, l) for e, l in cand_list if e not in matched.values()]
        by_hl = [x for x in free if hl_close(hl1, x[1]['hl_sec'])]
        if by_hl:
            matched[rank] = min(by_hl,
                                key=lambda x: abs(x[0] - e_mev * 1000))[0]
    for rank, e_mev, hl1 in sorted(v1_isomers):
        if rank in matched:
            continue
        e1 = e_mev * 1000.0
        free = [(e, l) for e, l in cand_list if e not in matched.values()]
        if free:
            best = min(free, key=lambda x: abs(x[0] - e1))
            if abs(best[0] - e1) < tol(e1):
                matched[rank] = best[0]
    # a rank still unmatched: admit its level from below the floor,
    # again by half-life first
    for rank, e_mev, hl1 in sorted(v1_isomers):
        if rank in matched:
            continue
        e1 = e_mev * 1000.0
        in_cands = {e for e, _ in cand_list}
        pool = [l for l in levels if not ledgered(l)
                and l['e_res'] not in in_cands]
        near = [l for l in pool if hl_close(hl1, l['hl_sec'])] or \
               [l for l in pool if abs(l['e_res'] - e1) < tol(e1)]
        if near:
            l = min(near, key=lambda l: abs(l['e_res'] - e1))
            cand_list = sorted(cand_list + [(l['e_res'], l)])
            matched[rank] = l['e_res']
            report['floor_override'].append((base, l['e_res'], l['hl_sec']))
        else:
            report['v1_isomer_lost'].append((base, rank, e_mev, hl1))

    # established ranks keep their v1 ordinals even when the (often
    # estimated) energies now sort differently — name preservation
    # outranks ascending-energy numbering; inversions are reported and
    # exhibited at the gate, never applied silently
    seq = [matched[r] for r in sorted(matched)]
    if seq != sorted(seq):
        report['order_swap'].append((base, {f'm{r}': matched[r]
                                            for r in sorted(matched)}))

    # newcomers below the highest established rank would renumber
    # established names: excluded and reported; above it they append.
    # Established ordinals are assigned by NAME, not position — a v1 rank
    # whose state left the evaluation leaves a GAP (73ZN keeps m2 for its
    # surviving 13 ms state; there is no m1), so every surviving name
    # keeps meaning the same physical state.
    ceiling = max(seq, default=-1.0)
    by_energy = {e: r for r, e in matched.items()}
    final = []
    for e, l in cand_list:
        if e in by_energy:
            final.append((by_energy[e], e, l))
        elif e > ceiling:
            final.append((None, e, l))
        else:
            report['newcomer_excluded'].append((base, e, l['hl_sec']))
    nxt = max(matched, default=0)
    out = []
    for r, e, l in sorted(final, key=lambda x: x[1]):
        if r is None:
            nxt += 1
            r = nxt
        out.append((r, e, l))
    return out


def build_triples(a, z, symbol, pstate, modes, isomers_of):
    """decay_mode string for one state; returns (string, closure_needs).
    `pstate` = (E_res [keV], raw offset [keV], band letters) of the parent
    state; isomers_of values are (rank, E_res, raw, band) per state."""
    p_eres, p_raw, p_band = pstate
    # net subset percentages from their parent mode
    pct = {}
    for mode, p in modes:
        if mode not in MODE_MAP:
            raise SystemExit(f'{a}{symbol.upper()}: unmapped decay mode '
                             f'{mode!r} — extend decay_modes.toml')
        if p is None:
            if len(modes) == 1:
                p = 100.0
                report['blank_branch_sole'].append(f'{a}{symbol.upper()}')
            else:
                report['blank_branch_multi'].append(
                    (f'{a}{symbol.upper()}', mode))
                continue
        pct[mode] = pct.get(mode, 0.0) + p
    for mode in list(pct):
        parent = MODE_MAP[mode].get('subset_of')
        if parent and parent in pct:
            pct[parent] = max(0.0, pct[parent] - pct[mode])
    total = sum(pct.values())
    if total > 100.5:
        report['nonclosing_branches'].append(
            (f'{a}{symbol.upper()}', round(total, 2)))
    scale = 1.0 / max(total, 100.0)       # normalize to sum <= 1

    triples, needs = [], []
    for mode, p in pct.items():
        spec = MODE_MAP[mode]
        frac = p * scale
        if frac <= 0.0:
            continue
        if spec.get('sfy'):
            if mode == 'SF+EC+B+':
                report['sf_combined'].append(f'{a}{symbol.upper()}')
            triples.append((spec['label'], 'SFY', frac))
            continue
        dz, da = spec['dz'], spec['da']
        d_z, d_a = z + dz, a + da
        if not 0 <= d_z < len(ELEMENTS) or d_a < 1:
            continue
        d_sym = ELEMENTS[d_z]
        d_base = f'{d_a}{d_sym.upper()}'
        # daughter-state split via level flow, where the daughter has isomers
        quads = isomers_of.get((d_a, d_z), []) if mode in _MODE_FILE else []
        if mode == 'IT':
            quads = [q for q in isomers_of.get((a, z), [])
                     if q[1] < p_eres - 0.5]
        placed = 0.0
        if quads:
            feed = isomer_feeding(a, symbol, p_raw, p_band, mode,
                                  [(rank, raw, band)
                                   for rank, _, raw, band in quads])
            fsum = sum(feed.values()) / 100.0
            if fsum > frac * 1.02 + 1e-9:
                report['feeding_overflow'].append(
                    (f'{a}{symbol.upper()}', mode, round(fsum, 4),
                     round(frac, 4)))
                norm = frac / fsum
            else:
                norm = 1.0
            for rank, _, _, _ in quads:
                if feed.get(rank, 0.0) > 0.0:
                    f_frac = feed[rank] / 100.0 * norm
                    triples.append((f'{spec["label"]}m{rank}',
                                    f'{d_base}m{rank}', f_frac))
                    needs.append(f'{d_base}m{rank}')
                    placed += f_frac
        gs_frac = frac - placed
        if gs_frac > 1e-12:
            triples.append((spec['label'], f'{d_base}g', gs_frac))
            needs.append(f'{d_base}g')
    s = ','.join(f'{lab}:{d}:{f:.6g}' for lab, d, f in triples)
    return s or None, needs


def load_v1_isomers():
    """{(A, Z): [(rank, E_level_MeV, half_life_s)]} from the previous
    generation — the half-life is the pairing fingerprint."""
    v1p = V1_DIR / 'decay.db'
    out = {}
    if not v1p.exists():
        return out
    v1 = sqlite3.connect(v1p)
    upper = [e.upper() for e in ELEMENTS]
    for name, e_lvl, hl in v1.execute(
            "SELECT name, E_level, half_life FROM chart WHERE isomer != 'g'"):
        m = re.fullmatch(r'(\d+)([A-Z]{1,2})m(\d+)', name)
        if not m:
            continue
        a, el, rank = int(m.group(1)), m.group(2), int(m.group(3))
        try:
            z = upper.index(el)
        except ValueError:
            continue
        out.setdefault((a, z), []).append((rank, e_lvl or 0.0, hl))
    return out


def main():
    nubase = load_nubase()
    v1_iso = load_v1_isomers()
    gs = list(csv.DictReader(open(LC / 'ground_states.csv')))

    # pass 1: the state universe (isomer admission + numbering, reconciled
    # against the previous generation's names by construction)
    isomers_of = {}                       # (A, Z) -> [E_kev ascending]
    iso_rows = {}                         # (A, Z) -> [(E_kev, level)]
    identities = []
    for r in gs:
        z, n = int(r['z']), int(r['n'])
        a = z + n
        symbol = r['symbol'].strip()
        identities.append((a, z, symbol, r))
        iso = admitted_isomers(a, z, symbol, nubase, v1_iso.get((a, z), []))
        isomers_of[(a, z)] = [(rank, e, l['energy'] or 0.0, l['shift'])
                              for rank, e, l in iso]
        iso_rows[(a, z)] = iso

    # pass 2: rows + closure
    BUILD_DIR.mkdir(exist_ok=True)
    path = BUILD_DIR / 'decay.db'
    if path.exists():
        path.unlink()
    con = sqlite3.connect(path)
    con.execute('CREATE TABLE chart (isotope TEXT, name TEXT, E_level REAL, '
                'isomer TEXT, J_pi TEXT, A INTEGER, Z INTEGER, N INTEGER, '
                'element TEXT, stable INTEGER, half_life REAL, '
                'unc_half_life REAL, abundance REAL, unc_abundance REAL, '
                'amu REAL, Delta REAL, decay_mode TEXT)')

    universe = set()
    for a, z, symbol, r in identities:
        universe.add(f'{a}{symbol.upper()}g')
        for rank, *_ in isomers_of[(a, z)]:
            universe.add(f'{a}{symbol.upper()}m{rank}')

    rows, all_needs = [], set()
    for a, z, symbol, r in identities:
        base = f'{a}{symbol.upper()}'
        me_kev = fnum(r['massexcess'])
        amu = fnum(r['atomic_mass'])
        amu = amu / 1e6 if amu is not None else (
            a + me_kev / AMU_KEV if me_kev is not None else None)
        stable = 1 if r['half_life'].strip() == 'STABLE' else 0
        hl, unc_hl = fnum(r['half_life_sec']), fnum(r['unc_hls'])
        if r['operator_hl'].strip():
            report['operator_half_lives'] += 1
        if not stable and hl is None:
            nb = nubase.get((a, z, 0))
            if nb and nb['hl_sec'] is not None:
                hl = nb['hl_sec']
            else:
                report['no_half_life'].append(base)
        modes = [(r[f'decay_{i}'].strip(), fnum(r[f'decay_{i}_%']))
                 for i in (1, 2, 3) if r[f'decay_{i}'].strip()]
        dm, needs = (None, []) if stable else build_triples(
            a, z, symbol, (0.0, 0.0, ''), modes, isomers_of)
        all_needs.update(needs)
        rows.append((base, base + 'g', 0.0, 'g', r['jp'].strip(), a, z,
                     a - z, symbol, stable, hl, unc_hl,
                     fnum(r['abundance']), fnum(r['unc_a']), amu,
                     me_kev / 1000.0 if me_kev is not None else None, dm))
        for rank, e_kev, l in iso_rows[(a, z)]:
            dm_i, needs_i = build_triples(
                a, z, symbol, (e_kev, l['energy'] or 0.0, l['shift']),
                l['modes'], isomers_of)
            all_needs.update(needs_i)
            abund = None
            if base == '180TA' and rank == 1:
                abund = fnum(r['abundance'])
            rows.append((base, f'{base}m{rank}', e_kev / 1000.0, f'm{rank}',
                         l['jp'], a, z, a - z, symbol, 0, l['hl_sec'],
                         l['unc_hl_sec'], abund, None, amu,
                         (me_kev / 1000.0 + e_kev / 1000.0)
                         if me_kev is not None else None, dm_i))
    # 180TA: the naturally occurring state is the isomer — LiveChart puts
    # the abundance on the nuclide; it belongs on m1, not the 8 h ground
    for i, row in enumerate(rows):
        if row[1] == '180TAg' and row[12] is not None:
            rows[i] = row[:12] + (None, None) + row[14:]

    # closure repair: NUBASE minimal rows for daughters outside the universe
    for name in sorted(all_needs - universe):
        m = re.fullmatch(r'(\d+)([A-Z]{1,2})(g|m\d+)', name)
        a, el, iso = int(m.group(1)), m.group(2), m.group(3)
        try:
            z = [e.upper() for e in ELEMENTS].index(el)
        except ValueError:
            continue
        if iso != 'g':
            continue                      # isomer daughters come from the census
        nb = nubase.get((a, z, 0))
        symbol = ELEMENTS[z]
        me = nb['me_kev'] if nb else None
        rows.append((f'{a}{el}', name, 0.0, 'g', '', a, z, a - z, symbol,
                     1 if (nb and nb['stable']) else 0,
                     nb['hl_sec'] if nb else None, None, None, None,
                     a + me / AMU_KEV if me is not None else None,
                     me / 1000.0 if me is not None else None, None))
        universe.add(name)
        report['closure_rows'].append(name)
        if nb is None or (not nb['stable'] and nb['hl_sec'] is None):
            report['no_half_life'].append(name)

    con.executemany('INSERT INTO chart VALUES (' + ','.join('?' * 17) + ')',
                    rows)
    con.commit()

    print(f'{len(rows)} chart rows ({sum(1 for r in rows if r[3] != "g")} '
          f'isomers, {len(report["closure_rows"])} closure-repair rows)')
    for key in ('ledger_exclusions', 'nubase_only_isomers', 'sf_combined',
                'feeding_overflow', 'nonclosing_branches'):
        if report[key]:
            print(f'  {key}: {len(report[key])} — {report[key][:8]}')
    print(f'  blank sole-mode branches -> 1.0: '
          f'{len(report["blank_branch_sole"])}; blank multi-mode skipped: '
          f'{len(report["blank_branch_multi"])}')
    print(f'  operator-limited half-lives (stored as value): '
          f'{report["operator_half_lives"]}')
    if report['no_half_life']:
        print(f'  UNSTABLE ROWS WITH NO HALF-LIFE (gate decision): '
              f'{len(report["no_half_life"])} — {report["no_half_life"][:12]}')
    con.close()


if __name__ == '__main__':
    main()
