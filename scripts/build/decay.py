"""Build decay.db v2 — the chart and radiation tables — from the LiveChart mirror.

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

Usage: python3 decay.py   # writes build/decay.db (chart + radiation)
"""
import csv
import math
import os
import re
import sqlite3
import sys
import tomllib
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import ARCHIVE_ROOT, BUILD_DIR, ELEMENTS, write_version

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
          'energy_estimated': [], 'order_swap': [],
          'blank_intensity_dropped': {}, 'orphan_state_rows': {},
          'xray_append_mismatch': [], 'unmapped_rad_decay': {},
          'norm_records_joined': 0}


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


def _jp_norm(s):
    return (s or '').replace('(', '').replace(')', '').replace(' ', '')


def match_state(row, states):
    """Resolve a radiation row to its parent state. `states` entries are
    (name, E_level [MeV], raw offset [keV], band, jp). The files key rows
    by p_energy WITHOUT the band letter (234Pa's X-band isomer keys as
    '73.92'), and two states can share an offset (110Rh keys both its
    ground state and its zero-offset band isomer as '0') — resolved by
    the row's J_pi."""
    p, p_band = level_token(row.get('p_energy'))
    if p is None:
        p, p_band = 0.0, ''
    hits = [s for s in states if abs(s[2] - p) < 1.0
            and (s[3] == p_band or not p_band)]
    if len(hits) > 1:
        rjp = _jp_norm(row.get('jp'))
        byjp = [s for s in hits if rjp and _jp_norm(s[4]) == rjp]
        if byjp:
            hits = byjp
        else:
            exact = [s for s in hits if s[3] == p_band]
            if exact:
                hits = exact
    return hits[0] if hits else None


# --- level-flow feeding: daughter-isomer split of one (state, mode) ---

_MODE_FILE = {'B-': ('decay_rads_bm.csv', ('intensity_beta',)),
              'EC+B+': ('decay_rads_bp.csv', ('intensity_beta', 'intensity_ec')),
              'EC': ('decay_rads_bp.csv', ('intensity_beta', 'intensity_ec')),
              'B+': ('decay_rads_bp.csv', ('intensity_beta', 'intensity_ec')),
              'A': ('decay_rads_a.csv', ('intensity',)),
              'IT': (None, ())}


def isomer_feeding(a, symbol, states, target, lc_mode, daughter_states):
    """Percent of `target`-state decays through `lc_mode` that populate
    each daughter isomer state: direct particle feeding plus gamma flow
    into the level weighted by (1 + conversion coefficient). Rows resolve
    to their parent state through match_state. Daughter states are
    (key, raw offset [keV], band letters) — a state in an unanchored band
    is fed by rows whose level token names the same band."""
    feed = {key: 0.0 for key, _, _ in daughter_states}

    def matched_keys(off, band):
        return [key for key, raw, b in daughter_states
                if b == band and abs(off - raw) < 1.0]

    def credit(token, amount):
        off, band = level_token(token)
        if off is None:
            return
        for key in matched_keys(off, band):
            feed[key] += amount

    fname, cols = _MODE_FILE[lc_mode]
    if fname:
        for r in rad_rows(a, symbol, fname):
            if match_state(r, states) is not target \
                    or r.get('decay', '').strip() != lc_mode:
                continue
            amount = sum(fnum(r.get(c)) or 0.0 for c in cols)
            credit(r.get('daughter_level_energy'), amount)
    for r in rad_rows(a, symbol, 'decay_rads_g.csv'):
        if match_state(r, states) is not target \
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
        if start is not None and matched_keys(start, s_band):
            continue
        cc = fnum(r.get('conversion_coeff')) or 0.0
        credit(r.get('end_level_energy'), i * (1.0 + cc))
    return feed


# --- chart assembly ---

def hl_close(h1, h2, factor=2.0):
    return (h1 and h2 and h1 > 0 and h2 > 0
            and abs(math.log10(h1 / h2)) < math.log10(factor))


def _level_energy(l, a, z, nubase, base):
    """Resolve a level's ordering/storage energy [keV]: LiveChart value,
    or — when the absolute energy is unknown (energy_shift set / zero for
    a non-ground level) — a NUBASE excitation estimate. The estimate is
    chosen by the level's half-life fingerprint first, energy closeness
    second: a banded level's raw offset can sit nearer another state's
    NUBASE entry than its own (203Pb's 122 ns level at 2923+X lies 26 keV
    from the 0.48 s state's 2949.2), and stealing that entry hands the
    established state's energy — and downstream its name and radiation
    rows — to the wrong level."""
    e = l['energy']
    if not l['shift'] and e:
        return e
    ents = [nubase[(a, z, i)] for i in (1, 2, 3, 4)
            if (a, z, i) in nubase and nubase[(a, z, i)]['exc_kev']]
    pool = [n for n in ents if hl_close(l['hl_sec'], n['hl_sec'])] or ents
    est = min(pool, key=lambda n: abs(n['exc_kev'] - (e or 0.0)))['exc_kev'] \
        if pool else None
    if est is not None:
        report['energy_estimated'].append((base, e, est))
        return est
    report.setdefault('energy_unresolved', []).append((base, e))
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

    cands = {}                            # keyed by LEVEL INDEX — two
    for l in levels:                      # levels can share a resolved
        if ledgered(l):                   # energy (85Nb's Y and Z bands)
            report['ledger_exclusions'].append((base, l['e_res']))
            continue
        if l['hl_sec'] is not None and l['hl_sec'] >= ISOMER_FLOOR_S:
            if any(abs(o['e_res'] - l['e_res']) < 1e-9
                   for o in cands.values()):
                report.setdefault('e_res_collision', []).append(
                    (base, l['e_res']))
            cands[l['idx']] = l
    i = 1
    while (a, z, i) in nubase:            # NUBASE census additions
        exc = nubase[(a, z, i)]['exc_kev']
        if exc is not None and not any(
                abs(exc - l['e_res']) < max(2.0, 0.01 * exc)
                for l in cands.values()):
            near = [l for l in levels
                    if abs(l['e_res'] - exc) < max(2.0, 0.01 * exc)
                    and not ledgered(l)]
            if near:
                l = min(near, key=lambda l: abs(l['e_res'] - exc))
                cands[l['idx']] = l
                report['nubase_only_isomers'].append((base, l['e_res']))
        i += 1
    cand_list = sorted(((l['e_res'], l) for l in cands.values()),
                       key=lambda x: (x[0], x[1]['idx']))
    if not v1_isomers:
        return [(i, e, l) for i, (e, l) in enumerate(cand_list, start=1)]

    # pair established ranks to candidates: the HALF-LIFE is the physical
    # fingerprint (the previous generation stored placeholder energies —
    # 0.1, 0.2 MeV — for states whose excitation was unknown, so energy
    # alone misleads); a factor-two half-life agreement matches, with
    # energy closeness as the tiebreak, and a plain energy match as the
    # fallback for ranks whose half-life was revised
    def tol(e):
        return max(30.0, 0.25 * e)

    # three passes in decreasing confidence: a tight energy agreement
    # outranks the half-life (the previous generation duplicated a state
    # under two ranks with one placeholder energy — 65FE — and the
    # duplicate must lose to the rank whose energy pins the level), the
    # half-life outranks a loose energy (placeholder energies are
    # arbitrary), and a loose energy catches revised half-lives.
    # Matches are held by LEVEL identity, not energy value — two
    # candidates can share a resolved energy (the 85Nb band pair), and a
    # value-keyed match would hand the rank to whichever sorts first and
    # knock both out of the later passes' free pool.
    matched = {}                          # v1 rank -> (energy, level)

    def free():
        used = {id(l) for _, l in matched.values()}
        return [(e, l) for e, l in cand_list if id(l) not in used]

    for rank, e_mev, hl1 in sorted(v1_isomers):
        e1 = e_mev * 1000.0
        tight = [x for x in free() if abs(x[0] - e1) < max(5.0, 0.01 * e1)]
        if tight:
            matched[rank] = min(tight, key=lambda x: abs(x[0] - e1))
    for rank, e_mev, hl1 in sorted(v1_isomers):
        if rank in matched:
            continue
        by_hl = [x for x in free() if hl_close(hl1, x[1]['hl_sec'])]
        if by_hl:
            matched[rank] = min(by_hl,
                                key=lambda x: abs(x[0] - e_mev * 1000))
    for rank, e_mev, hl1 in sorted(v1_isomers):
        if rank in matched:
            continue
        e1 = e_mev * 1000.0
        pool = free()
        if pool:
            best = min(pool, key=lambda x: abs(x[0] - e1))
            if abs(best[0] - e1) < tol(e1):
                matched[rank] = best
    # a rank still unmatched: admit its level from below the floor,
    # again by half-life first
    for rank, e_mev, hl1 in sorted(v1_isomers):
        if rank in matched:
            continue
        e1 = e_mev * 1000.0
        in_cands = {id(l) for _, l in cand_list}
        pool = [l for l in levels if not ledgered(l)
                and id(l) not in in_cands]
        near = [l for l in pool if hl_close(hl1, l['hl_sec'])] or \
               [l for l in pool if abs(l['e_res'] - e1) < tol(e1)]
        if near:
            l = min(near, key=lambda l: abs(l['e_res'] - e1))
            cand_list = sorted(cand_list + [(l['e_res'], l)],
                               key=lambda x: (x[0], x[1]['idx']))
            matched[rank] = (l['e_res'], l)
            report['floor_override'].append((base, l['e_res'], l['hl_sec']))
        else:
            report['v1_isomer_lost'].append((base, rank, e_mev, hl1))

    # established ranks keep their v1 ordinals even when the (often
    # estimated) energies now sort differently — name preservation
    # outranks ascending-energy numbering; inversions are reported and
    # exhibited at the gate, never applied silently
    seq = [matched[r][0] for r in sorted(matched)]
    if seq != sorted(seq):
        report['order_swap'].append((base, {f'm{r}': matched[r][0]
                                            for r in sorted(matched)}))

    # newcomers below the highest established rank would renumber
    # established names: excluded and reported; above it they append.
    # Established ordinals are assigned by NAME, not position — a v1 rank
    # whose state left the evaluation leaves a GAP (73ZN keeps m2 for its
    # surviving 13 ms state; there is no m1), so every surviving name
    # keeps meaning the same physical state.
    ceiling = max(seq, default=-1.0)
    rank_of = {id(l): r for r, (e, l) in matched.items()}
    final = []
    for e, l in cand_list:
        if id(l) in rank_of:
            final.append((rank_of.pop(id(l)), e, l))
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


def states_of(gs_row, iso):
    """The state list for one nuclide: (name, E_level [MeV], raw offset
    [keV], band, jp) — ground state first."""
    return [('g', 0.0, 0.0, '', gs_row['jp'].strip())] + [
        (f'm{rank}', e_res / 1000.0, l['energy'] or 0.0, l['shift'], l['jp'])
        for rank, e_res, l in iso]


def build_triples(a, z, symbol, states, target, modes, isomers_of):
    """decay_mode string for one state; returns (string, closure_needs).
    `target` is the parent state's entry in `states`; isomers_of values
    are (rank, E_res, raw, band) per state."""
    p_eres = target[1] * 1000.0
    # net subset percentages from their parent mode. An IT entry on a
    # GROUND state is a bookkeeping artifact of historically swapped
    # state assignments (60Mn) — dropped, with the remaining modes
    # renormalized to carry the full decay.
    gs_it = [p for m, p in modes if m == 'IT'] if target[0] == 'g' else []
    if gs_it:
        report.setdefault('gs_it_dropped', []).append(f'{a}{symbol.upper()}')
        modes = [(m, p) for m, p in modes if m != 'IT']
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
        if not parent:
            continue
        candidates = [parent] + (['EC', 'B+'] if parent == 'EC+B+' else [])
        for cand in candidates:
            if cand in pct:
                pct[cand] = max(0.0, pct[cand] - pct[mode])
                break
    total = sum(pct.values())
    if total > 100.5:
        report['nonclosing_branches'].append(
            (f'{a}{symbol.upper()}', round(total, 2)))
    # normalize to sum <= 1; when a ground-state IT artifact was dropped,
    # the surviving modes carry the whole decay
    scale = 1.0 / total if (gs_it and total > 0) else 1.0 / max(total, 100.0)

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
            report.setdefault('off_chart_daughter', []).append(
                (f'{a}{symbol.upper()}', mode))
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
            feed = isomer_feeding(a, symbol, states, target, mode,
                                  [(rank, raw, band)
                                   for rank, _, raw, band in quads])
            fsum = sum(feed.values()) / 100.0
            if fsum > frac * 1.02 + 1e-9:
                report['feeding_overflow'].append(
                    (f'{a}{symbol.upper()}', mode, round(fsum, 4),
                     round(frac, 4)))
            # the split can never exceed the mode total (round-off and
            # relative-intensity datasets would otherwise push a branch
            # over 1)
            norm = min(1.0, frac / fsum) if fsum > 0 else 1.0
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
    s_tot = sum(f for _, _, f in triples)
    if s_tot > 1.0:                       # round-off from the split
        triples = [(lab, d, f / s_tot) for lab, d, f in triples]
    s = ','.join(f'{lab}:{d}:{f:.6g}' for lab, d, f in triples)
    return s or None, needs


# --- radiation tables ---

_XR_NOTES = {'KA1': 'XR k_alpha_1', 'KA2': 'XR k_alpha_2',
             'KpB1': 'XR k_beta_1', 'KpB2': 'XR k_beta_2', 'L': 'XR l'}
_CE_SHELLS = ('K', 'L', 'M', 'N', 'N+', 'O', 'P')


def _rad_label(lc_decay):
    """LiveChart 'decay' string -> stored decay_mode label."""
    spec = MODE_MAP.get(lc_decay)
    if spec:
        return spec['label']
    report['unmapped_rad_decay'][lc_decay] = \
        report['unmapped_rad_decay'].get(lc_decay, 0) + 1
    return lc_decay


def _drop_blank(table, base):
    counts = report['blank_intensity_dropped'].setdefault(table, {})
    counts[base] = counts.get(base, 0) + 1



def _accept_intensity(i, table, base, energy=None):
    """Shared intensity gate: blank/nonpositive rows drop with the blank
    accounting; rows above 1000 percent are relative-intensity outliers
    and drop with their own accounting."""
    if i is None or i <= 0:
        _drop_blank(table, base)
        return False
    if i > 1000.0:
        report.setdefault('over_1000_dropped', []).append(
            (base, table, energy, i))
        return False
    return True


def _gamma_row(base, st, energy, i, unc, note, lab, norms):
    return (base, st[1], st[0], energy, i, unc, note, lab,
            _state_norm(norms, base, st, lab))


def _state_norm(norms, base, st, lab):
    """unc_norm join. A zero-offset band isomer would collide with the
    ground state's (base, 0.0, mode) key — no join for that class."""
    if st[3] and st[2] == 0.0:
        return None
    return norms.get((base, round(st[2], 1), lab))


def build_radiation(con, identities, iso_rows, norms):
    """The five radiation tables, keyed (isotope short name, isomer) with
    intensities as served: absolute percent per 100 decays of the parent
    state, branching folded in. gammas carry the additive unc_norm column
    (relative NR*BR normalization uncertainty of the line's dataset)."""
    con.execute('CREATE TABLE gammas (isotope TEXT, E_level REAL, '
                'isomer TEXT, energy REAL, intensity REAL, '
                'unc_intensity REAL, notes TEXT, decay_mode TEXT, '
                'unc_norm REAL)')
    con.execute('CREATE TABLE electrons (isotope TEXT, E_level REAL, '
                'isomer TEXT, energy REAL, intensity REAL, '
                'unc_intensity REAL, notes TEXT, decay_mode TEXT)')
    con.execute('CREATE TABLE alphas (isotope TEXT, E_level REAL, '
                'isomer TEXT, energy REAL, intensity REAL, '
                'unc_intensity REAL, decay_mode TEXT)')
    for t in ('beta_minus', 'beta_plus'):
        con.execute(f'CREATE TABLE {t} (isotope TEXT, E_level REAL, '
                    f'isomer TEXT, mean_energy REAL, intensity REAL, '
                    f'unc_intensity REAL, endpoint_energy REAL, '
                    f'decay_mode TEXT)')
    con.execute('CREATE TABLE elements (Z INTEGER, element TEXT)')
    con.executemany('INSERT INTO elements VALUES (?,?)',
                    [(z, ELEMENTS[z]) for z in range(1, len(ELEMENTS))])

    g_rows, e_rows, a_rows, bm_rows, bp_rows = [], [], [], [], []
    for a, z, symbol, r in identities:
        base = f'{a}{symbol.upper()}'
        states = states_of(r, iso_rows[(a, z)])

        def find_state(row):
            st = match_state(row, states)
            if st is None:
                report['orphan_state_rows'][base] = \
                    report['orphan_state_rows'].get(base, 0) + 1
            return st

        # gammas: the g file's PLACED rows only. Its unplaced rows are the
        # x-ray append (cross-checked against the x file below) plus the
        # 511 keV annihilation rows, which are taken from the file as
        # served (synthesis from the beta-plus total is the fallback).
        unplaced = []
        file_511 = set()
        for row in rad_rows(a, symbol, 'decay_rads_g.csv'):
            placed = (row.get('start_level_energy', '').strip()
                      or row.get('end_level_energy', '').strip())
            e_g = fnum(row.get('energy'))
            i = fnum(row.get('intensity'))
            if not placed and e_g is not None and abs(e_g - 511.0) < 0.5:
                st = find_state(row)
                if st is None or i is None:
                    continue
                lab = _rad_label(row.get('decay', '').strip())
                # only positron-capable datasets produce annihilation: a
                # beta-minus line near 511 keV (142Cs) keeps its own
                # energy as an unplaced gamma instead
                if lab != 'E':
                    unplaced.append(row)
                    continue
                g_rows.append(_gamma_row(base, st, 511.0, i,
                                         fnum(row.get('unc_i')),
                                         'Annihil.', lab, norms))
                file_511.add((st, lab))
                continue
            if not placed:
                unplaced.append(row)
                continue
            st = find_state(row)
            if st is None:
                continue
            if not _accept_intensity(i, 'gammas', base, e_g):
                continue
            lab = _rad_label(row.get('decay', '').strip())
            g_rows.append(_gamma_row(base, st, e_g, i,
                                     fnum(row.get('unc_i')), '', lab, norms))

        # x-rays: exclusively from the x file; KB duplicates KpB1 and the
        # L-subshell lines duplicate the aggregate L row where it exists
        by_state_mode = {}
        x_all = []
        for row in rad_rows(a, symbol, 'decay_rads_x.csv'):
            x_all.append((fnum(row.get('energy')), fnum(row.get('intensity'))))
            st = find_state(row)
            if st is None:
                continue
            key = (st, row.get('decay', '').strip())
            by_state_mode.setdefault(key, []).append(row)
        for (st, dmode), rows in by_state_mode.items():
            shells = {r.get('shell', '').strip() for r in rows}
            lab = _rad_label(dmode)
            for row in rows:
                sh = row.get('shell', '').strip()
                if sh == 'KB':
                    continue              # byte-duplicate of KpB1
                if sh.startswith('L') and len(sh) > 1 and 'L' in shells:
                    continue              # subshell under an aggregate L
                note = _XR_NOTES.get(sh)
                if note is None:
                    if sh.startswith('L'):
                        note = 'XR l'     # sub-resolved, no aggregate
                    else:
                        raise SystemExit(f'{base}: unmapped x-ray shell '
                                         f'{sh!r}')
                i = fnum(row.get('intensity'))
                if not _accept_intensity(i, 'gammas', base,
                                         fnum(row.get('energy'))):
                    continue
                g_rows.append(_gamma_row(base, st, fnum(row.get('energy')),
                                         i, fnum(row.get('unc_i')), note,
                                         lab, norms))
        # the g file's unplaced rows are the x-ray append plus genuinely
        # UNPLACED gammas (observed lines without a level assignment) —
        # consume the x-file rows from the unplaced multiset; whatever
        # remains is a real gamma and is kept; an x-file row the append
        # never contained is a source inconsistency, reported
        x_left = Counter((e, i) for e, i in x_all if i is not None)
        for row in unplaced:
            e_g = fnum(row.get('energy'))
            i = fnum(row.get('intensity'))
            if i is not None and x_left.get((e_g, i), 0) > 0:
                x_left[(e_g, i)] -= 1
                continue                  # the x-ray append copy
            st = find_state(row)
            if st is None:
                continue
            if not _accept_intensity(i, 'gammas', base, e_g):
                continue
            report['unplaced_gammas_kept'] = \
                report.get('unplaced_gammas_kept', 0) + 1
            lab = _rad_label(row.get('decay', '').strip())
            g_rows.append(_gamma_row(base, st, e_g, i,
                                     fnum(row.get('unc_i')), '', lab, norms))
        if +x_left:
            report['xray_append_mismatch'].append(base)

        # electrons: conversion electrons + Auger; the AU K aggregate is
        # the exact sum of its KLL/KLX/KXY subshells, which then drop
        e_groups = {}
        for row in rad_rows(a, symbol, 'decay_rads_e.csv'):
            st = find_state(row)
            if st is None:
                continue
            e_groups.setdefault((st, row.get('decay', '').strip()),
                                []).append(row)
        for (st, dmode), rows in e_groups.items():
            lab = _rad_label(dmode)
            au_shells = {r.get('shell', '').strip() for r in rows
                         if r.get('type', '').strip() == 'AU'}
            for row in rows:
                typ = row.get('type', '').strip()
                sh = row.get('shell', '').strip()
                if typ == 'CE':
                    if sh not in _CE_SHELLS:
                        raise SystemExit(f'{base}: unmapped CE shell {sh!r}')
                    note = f'CE {"N" if sh == "N+" else sh}'
                elif typ == 'AU':
                    if sh in ('KLL', 'KLX', 'KXY') and 'K' in au_shells:
                        continue
                    note = f'Auger {sh[0]}'
                else:
                    raise SystemExit(f'{base}: unmapped electron type {typ!r}')
                i = fnum(row.get('intensity'))
                if not _accept_intensity(i, 'electrons', base,
                                         fnum(row.get('energy'))):
                    continue
                e_rows.append((base, st[1], st[0], fnum(row.get('energy')),
                               i, fnum(row.get('unc_i')), note, lab))

        # alphas
        for row in rad_rows(a, symbol, 'decay_rads_a.csv'):
            st = find_state(row)
            if st is None:
                continue
            i = fnum(row.get('intensity'))
            if not _accept_intensity(i, 'alphas', base,
                                     fnum(row.get('energy'))):
                continue
            a_rows.append((base, st[1], st[0], fnum(row.get('energy')), i,
                           fnum(row.get('unc_i')),
                           _rad_label(row.get('decay', '').strip())))

        # betas: bm rows as beta-minus; bp rows carry the beta-plus
        # component only (a blank intensity_beta is an EC-only transition,
        # not a beta_plus row). The 511 keV annihilation line is
        # synthesized per (state, mode) at twice the total beta-plus
        # intensity, as the previous generation stored it.
        for row in rad_rows(a, symbol, 'decay_rads_bm.csv'):
            st = find_state(row)
            if st is None:
                continue
            i = fnum(row.get('intensity_beta'))
            if not _accept_intensity(i, 'beta_minus', base,
                                     fnum(row.get('mean_energy'))):
                continue
            bm_rows.append((base, st[1], st[0], fnum(row.get('mean_energy')),
                            i, fnum(row.get('unc_ib')),
                            fnum(row.get('max_energy')),
                            _rad_label(row.get('decay', '').strip())))
        bp_tot = {}
        for row in rad_rows(a, symbol, 'decay_rads_bp.csv'):
            st = find_state(row)
            if st is None:
                continue
            i = fnum(row.get('intensity_beta'))
            if i is None:
                continue                  # EC-only transition
            if not _accept_intensity(i, 'beta_plus', base,
                                     fnum(row.get('mean_energy'))):
                continue
            lab = _rad_label(row.get('decay', '').strip())
            e_ec = fnum(row.get('energy_ec'))
            # the bp file serves no max_energy; the endpoint is the
            # transition Q_EC minus twice the electron mass
            endpoint = e_ec - 1022.0 if e_ec is not None else None
            bp_rows.append((base, st[1], st[0], fnum(row.get('mean_energy')),
                            i, fnum(row.get('unc_ib')), endpoint, lab))
            tot, var = bp_tot.get((st, lab), (0.0, 0.0))
            bp_tot[(st, lab)] = (tot + i,
                                 var + (fnum(row.get('unc_ib')) or 0.0) ** 2)
        for (st, lab), (tot, var) in bp_tot.items():
            if tot > 0 and (st, lab) not in file_511:
                report['annihil_synthesized'] = \
                    report.get('annihil_synthesized', 0) + 1
                g_rows.append(_gamma_row(base, st, 511.0, 2.0 * tot,
                                         2.0 * var ** 0.5 or None,
                                         'Annihil.', lab, norms))

    con.executemany('INSERT INTO gammas VALUES (?,?,?,?,?,?,?,?,?)', g_rows)
    con.executemany('INSERT INTO electrons VALUES (?,?,?,?,?,?,?,?)', e_rows)
    con.executemany('INSERT INTO alphas VALUES (?,?,?,?,?,?,?)', a_rows)
    con.executemany('INSERT INTO beta_minus VALUES (?,?,?,?,?,?,?,?)', bm_rows)
    con.executemany('INSERT INTO beta_plus VALUES (?,?,?,?,?,?,?,?)', bp_rows)
    con.commit()
    print(f'radiation: {len(g_rows)} gammas ({sum(1 for r in g_rows if r[8] is not None)} '
          f'with unc_norm), {len(e_rows)} electrons, {len(a_rows)} alphas, '
          f'{len(bm_rows)} beta-, {len(bp_rows)} beta+')


def load_norms(identities, iso_rows):
    """{(base, state_E_kev rounded, label): relative NR*BR uncertainty}
    from the ENSDF normalization records via paceENSDF."""
    try:
        import paceENSDF
    except ImportError:
        print('  (paceENSDF unavailable — unc_norm ships NULL)')
        return {}
    import contextlib
    import io
    e = paceENSDF.ENSDF()
    with contextlib.redirect_stdout(io.StringIO()):
        data = e.load_ensdf()
    norms = {}
    for a, z, symbol, r in identities:
        base = f'{a}{symbol.upper()}'
        pid = f'{symbol}{a}'
        n_states = 1 + len(iso_rows[(a, z)])
        # one index beyond our state count: paceENSDF's decay indices
        # count ITS datasets, which can exceed our admitted states (an
        # excluded level with its own dataset); a missing index returns
        # no record and costs nothing
        for idx in range(n_states + 1):
            for kw, lab in (('BM', 'B'), ('ECBP', 'E'), ('A', 'A')):
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        rec = e.norm_record(data, pid, idx, mode=kw)
                except Exception:
                    rec = None
                if not rec:
                    continue
                for key, vals in rec.items():
                    try:                  # band-level parents ('0.0+X') skip
                        e_kev = float(key[4])
                    except (TypeError, ValueError):
                        continue
                    nr, u_nr, ok_nr = vals[0]
                    br, u_br, ok_br = vals[2]
                    if not ok_nr or nr <= 0:
                        continue
                    rel = ((u_nr / nr) ** 2 +
                           ((u_br / br) ** 2 if br > 0 else 0.0)) ** 0.5
                    norms[(base, round(e_kev, 1), lab)] = rel
                    report['norm_records_joined'] += 1
    return norms


def nubase_iso_hl(nubase, a, z, e_kev):
    """NUBASE half-life for the isomer nearest `e_kev`, used when the
    LiveChart level carries none."""
    i = 1
    while (a, z, i) in nubase:
        nb = nubase[(a, z, i)]
        if (nb['exc_kev'] is not None and nb['hl_sec'] is not None
                and abs(nb['exc_kev'] - e_kev) < max(5.0, 0.02 * e_kev)):
            return nb['hl_sec']
        i += 1
    return None


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
            a + me_kev / AMU_KEV if me_kev is not None else float(a))
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
        states = states_of(r, iso_rows[(a, z)])
        dm, needs = (None, []) if stable else build_triples(
            a, z, symbol, states, states[0], modes, isomers_of)
        all_needs.update(needs)
        rows.append((base, base + 'g', 0.0, 'g', r['jp'].strip(), a, z,
                     a - z, symbol, stable, hl, unc_hl,
                     fnum(r['abundance']), fnum(r['unc_a']), amu,
                     me_kev / 1000.0 if me_kev is not None else None, dm))
        for i, (rank, e_kev, l) in enumerate(iso_rows[(a, z)]):
            dm_i, needs_i = build_triples(
                a, z, symbol, states, states[1 + i], l['modes'], isomers_of)
            all_needs.update(needs_i)
            abund = None
            if base == '180TA' and rank == 1:
                abund = fnum(r['abundance'])
            hl_i = l['hl_sec']
            if l['operator_hl']:
                report['operator_half_lives'] += 1
            if hl_i is None:              # census/override levels may lack one
                hl_i = nubase_iso_hl(nubase, a, z, e_kev)
                if hl_i is not None:
                    report.setdefault('isomer_hl_from_nubase', []).append(
                        f'{base}m{rank}')
            if hl_i is None:              # the state's own dataset rows
                for fname in ('decay_rads_g.csv', 'decay_rads_a.csv',
                              'decay_rads_bm.csv', 'decay_rads_bp.csv'):
                    for row in rad_rows(a, symbol, fname):
                        if match_state(row, states) is states[1 + i]:
                            hl_i = fnum(row.get('half_life_sec'))
                            if hl_i is not None:
                                break
                    if hl_i is not None:
                        report.setdefault('isomer_hl_from_dataset',
                                          []).append(f'{base}m{rank}')
                        break
            if hl_i is None:
                report['no_half_life'].append(f'{base}m{rank}')
            rows.append((base, f'{base}m{rank}', e_kev / 1000.0, f'm{rank}',
                         l['jp'], a, z, a - z, symbol, 0, hl_i,
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
            raise SystemExit(f'closure needs isomer {name} outside the '
                             f'universe — the split should only name '
                             f'admitted isomers')
        nb = nubase.get((a, z, 0))
        symbol = ELEMENTS[z]
        me = nb['me_kev'] if nb else None
        rows.append((f'{a}{el}', name, 0.0, 'g', '', a, z, a - z, symbol,
                     1 if (nb and nb['stable']) else 0,
                     nb['hl_sec'] if nb else None, None, None, None,
                     a + me / AMU_KEV if me is not None else float(a),
                     me / 1000.0 if me is not None else None, None))
        universe.add(name)
        report['closure_rows'].append(name)
        if nb is None or (not nb['stable'] and nb['hl_sec'] is None):
            report['no_half_life'].append(name)

    con.executemany('INSERT INTO chart VALUES (' + ','.join('?' * 17) + ')',
                    rows)
    con.commit()

    norms = load_norms(identities, iso_rows)
    build_radiation(con, identities, iso_rows, norms)
    write_version(con, 'decay',
                  'ENSDF via IAEA LiveChart API mirror (chart + radiation), '
                  'NUBASE2020/AME2020 closure and estimates, ENSDF N-records '
                  'via paceENSDF (unc_norm)',
                  'LiveChart mirror 2026-07-16, per-row ENSDF cut-offs '
                  'recorded in the source archive')
    con.commit()

    print(f'{len(rows)} chart rows ({sum(1 for r in rows if r[3] != "g")} '
          f'isomers, {len(report["closure_rows"])} closure-repair rows)')
    for key in ('ledger_exclusions', 'nubase_only_isomers', 'sf_combined',
                'feeding_overflow', 'nonclosing_branches',
                'newcomer_excluded', 'floor_override', 'v1_isomer_lost',
                'order_swap', 'gs_it_dropped', 'isomer_hl_from_nubase',
                'isomer_hl_from_dataset', 'over_1000_dropped'):
        if report.get(key):
            print(f'  {key}: {len(report[key])} — {report[key][:8]}')
    if report['energy_estimated']:
        print(f'  band/unknown level energies estimated from NUBASE: '
              f'{len(report["energy_estimated"])}; no estimate available '
              f'(raw offset or 0 kept): '
              f'{len(report.get("energy_unresolved", []))}')
    print(f'  blank sole-mode branches -> 1.0: '
          f'{len(report["blank_branch_sole"])}; blank multi-mode skipped: '
          f'{len(report["blank_branch_multi"])}')
    print(f'  operator-limited half-lives (stored as value): '
          f'{report["operator_half_lives"]}')
    if report['no_half_life']:
        print(f'  UNSTABLE ROWS WITH NO HALF-LIFE (gate decision): '
              f'{len(report["no_half_life"])} — {report["no_half_life"][:12]}')
    for table, counts in sorted(report['blank_intensity_dropped'].items()):
        print(f'  blank-intensity rows dropped from {table}: '
              f'{sum(counts.values())} across {len(counts)} nuclides')
    if report['orphan_state_rows']:
        print(f'  rows keyed to non-admitted parent states (dropped): '
              f'{sum(report["orphan_state_rows"].values())} across '
              f'{len(report["orphan_state_rows"])} nuclides — '
              f'{sorted(report["orphan_state_rows"], key=report["orphan_state_rows"].get)[-6:]}')
    if report['xray_append_mismatch']:
        print(f'  g-file x-ray append != x file: '
              f'{len(report["xray_append_mismatch"])} nuclides — '
              f'{report["xray_append_mismatch"][:8]}')
    if report['unmapped_rad_decay']:
        print(f'  radiation rows with unmapped decay strings (label kept '
              f'verbatim): {report["unmapped_rad_decay"]}')
    print(f'  unplaced gammas kept: {report.get("unplaced_gammas_kept", 0)}; '
          f'511 keV rows synthesized (no file row): '
          f'{report.get("annihil_synthesized", 0)}')
    print(f'  ENSDF normalization records joined: '
          f'{report["norm_records_joined"]}')
    con.close()


if __name__ == '__main__':
    main()
