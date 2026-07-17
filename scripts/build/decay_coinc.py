"""Add the coincidences table to decay.db.

Per parent state: the probability per 100 decays that a pair of photons
is emitted together — the input for true-coincidence-summing corrections
(summing-out via the rows keyed on the line being corrected, summing-in
via a self-join on pair energies). Source: the paceENSDF gamma-gamma and
gamma-X-ray coincidence structures (conversion, level feeding and decay
branching already folded by construction), plus 511 keV annihilation
pairs synthesized from the state's stored annihilation intensity.

    coincidences(isotope TEXT, isomer TEXT, energy REAL,
                 coinc_energy REAL, coinc_type TEXT,
                 intensity REAL, unc_intensity REAL)

Conventions: keying matches the radiation tables; `energy` values take
the byte-identical gammas-table energy of the same state via nearest-row
matching within 0.5 keV, so downstream SQL joins are exact (a pace line
with no gammas row in reach keeps its own energy and is counted);
gamma-gamma pairs store BOTH directions (mirrors added only where pace
did not itself serve the reverse direction; same-key entries sum); gamma-X and
gamma-511 pairs store both directions too, with the X or 511 partner
identified by coinc_type ('XR', '511'; plain gammas 'g'). The 511 pairs
use the same two-photons-per-beta-plus convention as the radiation
tables: p(gamma & 511) = I_gamma x I_511 / (100 x BR_ECBP); the
(511, 511) self-pair carries I_511/2 (one annihilation pair per beta
plus). Rows below 1e-4 percent are dropped (build parameter).

Usage: python3 decay_coinc.py       # appends to build/decay.db
"""
import contextlib
import io
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import BUILD_DIR

THRESHOLD_PCT = 1e-4
SNAP_KEV = 0.5
MODE_LABEL = {'betaMinusDecay': 'B', 'electronCaptureBetaPlusDecay': 'E',
              'alphaDecay': 'A'}

co_report = {'datasets_matched': 0, 'datasets_unmatched': [],
             'snap_missed': 0, 'threshold_dropped': 0,
             'over_singles': 0}


def chart_states(con):
    """{(Z, A): [(isomer, name-short, E_level_kev)]}"""
    out = {}
    for iso_short, isomer, a, z, e_lvl in con.execute(
            'SELECT isotope, isomer, A, Z, E_level FROM chart'):
        out.setdefault((z, a), []).append((isomer, iso_short,
                                           (e_lvl or 0.0) * 1000.0))
    return out


def gamma_index(con):
    """{(isotope, isomer, label): [(energy, intensity, notes)]} for the
    energy snap and the 511/BR lookups."""
    out = {}
    for iso, isomer, e, i, notes, lab in con.execute(
            'SELECT isotope, isomer, energy, intensity, notes, decay_mode '
            'FROM gammas'):
        out.setdefault((iso, isomer, lab), []).append((e, i, notes or ''))
    return out


def snap(gidx, iso, isomer, lab, e_kev, want_xr=False):
    rows = [r for r in gidx.get((iso, isomer, lab), [])
            if ('XR' in r[2]) == want_xr and r[2] != 'Annihil.']
    if rows:
        best = min(rows, key=lambda r: abs(r[0] - e_kev))
        if abs(best[0] - e_kev) < SNAP_KEV:
            return best[0], best[1]
    co_report['snap_missed'] += 1
    return e_kev, None


def ecbp_branch(con, iso, isomer):
    """The EC+beta+ branching fraction of the state, from its triples."""
    row = con.execute('SELECT decay_mode FROM chart WHERE isotope=? AND '
                      'isomer=?', (iso, isomer)).fetchone()
    if not row or not row[0]:
        return 0.0
    tot = 0.0
    for entry in row[0].split(','):
        parts = entry.split(':')
        if len(parts) == 3 and (parts[0] == 'E'
                                or parts[0].startswith('Em')):
            tot += float(parts[2])
    return tot


def main():
    import paceENSDF
    con = sqlite3.connect(BUILD_DIR / 'decay.db')
    if con.execute("SELECT name FROM sqlite_master WHERE name='coincidences'"
                   ).fetchone():
        con.execute('DROP TABLE coincidences')
    con.execute('CREATE TABLE coincidences (isotope TEXT, isomer TEXT, '
                'energy REAL, coinc_energy REAL, coinc_type TEXT, '
                'intensity REAL, unc_intensity REAL)')
    states = chart_states(con)
    gidx = gamma_index(con)

    e = paceENSDF.ENSDF()
    with contextlib.redirect_stdout(io.StringIO()):
        data = e.load_pace()

    rows = {}                             # (iso, isomer, e1, e2, typ) ->
    mirror = {}                           #   [intensity, var]; mirrors held
                                          #   back until the pass ends

    def add(iso, isomer, e1, e2, typ, p, up, mirrored=False):
        if p is None:
            co_report['missing_intensity'] = \
                co_report.get('missing_intensity', 0) + 1
            return
        if p < THRESHOLD_PCT:
            co_report['threshold_dropped'] += 1
            return
        key = (iso, isomer, e1, e2, typ)
        target = mirror if mirrored else rows
        if key in target:
            co_report['aggregated_pairs'] = \
                co_report.get('aggregated_pairs', 0) + 1
        acc = target.setdefault(key, [0.0, 0.0])
        acc[0] += p
        acc[1] += (up or 0.0) ** 2

    for d in data:
        lab = MODE_LABEL[d['decayMode']]
        cands = states.get((d['parentZ'], d['parentA']), [])
        try:
            e_par = float(d['parentDecayLevelEnergy'] or 0.0)
        except (TypeError, ValueError):   # band-level parents ('0.0+X')
            co_report.setdefault('datasets_band_parent', []).append(
                (d['parentID'], d['decayIndex'],
                 d['parentDecayLevelEnergy']))
            continue
        near = [c for c in cands if abs(c[2] - e_par) < 1.0]
        if not near:
            co_report['datasets_unmatched'].append(
                (d['parentID'], d['decayIndex'], e_par))
            continue
        isomer, iso, _ = near[0]
        co_report['datasets_matched'] += 1

        if d['datasetID'] == 'GG':
            for c in d.get('decayCoincidences', []):
                p = c['absoluteCoincidenceIntensity']
                up = c['dAbsoluteCoincidenceIntensity']
                e1, i1 = snap(gidx, iso, isomer, lab, c['gammaEnergyGate'])
                e2, i2 = snap(gidx, iso, isomer, lab,
                              c['gammaEnergyCoincidence'])
                # a pair intensity cannot exceed EITHER partner's singles
                if (i1 is not None and p > 1.02 * i1 + 0.01) or \
                        (i2 is not None and p > 1.02 * i2 + 0.01):
                    co_report['over_singles'] += 1
                add(iso, isomer, e1, e2, 'g', p, up)
                add(iso, isomer, e2, e1, 'g', p, up, mirrored=True)
        else:                             # GX
            for c in d.get('gammaXrayCoincidences', []):
                p = c['absoluteCoincidenceIntensityGammaXray']
                up = c['dAbsoluteCoincidenceIntensityGammaXray']
                e1, _ = snap(gidx, iso, isomer, lab, c['gammaEnergy'])
                ex, _ = snap(gidx, iso, isomer, lab, c['XrayEnergy'],
                             want_xr=True)
                add(iso, isomer, e1, ex, 'XR', p, up)
                add(iso, isomer, ex, e1, 'XR', p, up, mirrored=True)

    # 511 pairs from the stored annihilation intensities
    for (iso, isomer, lab), grows in gidx.items():
        if lab != 'E':
            continue
        i511 = sum(i for _, i, notes in grows if notes == 'Annihil.')
        if i511 <= 0:
            continue
        br = ecbp_branch(con, iso, isomer)
        if br <= 0:
            continue
        for e_g, i_g, notes in grows:
            if notes == 'Annihil.' or 'XR' in notes:
                continue
            p = i_g * i511 / (100.0 * br)
            add(iso, isomer, e_g, 511.0, '511', p, None)
            add(iso, isomer, 511.0, e_g, '511', p, None, mirrored=True)
        add(iso, isomer, 511.0, 511.0, '511', i511 / 2.0, None)

    for key, acc in mirror.items():
        if key not in rows:
            rows[key] = acc
    con.executemany(
        'INSERT INTO coincidences VALUES (?,?,?,?,?,?,?)',
        sorted((k[0], k[1], k[2], k[3], k[4], v[0],
                v[1] ** 0.5 if v[1] > 0 else None)
               for k, v in rows.items()))
    # the table is ~20x larger than the radiation tables; the summing
    # queries key on (isotope, isomer, energy)
    con.execute('CREATE INDEX idx_coinc ON coincidences '
                '(isotope, isomer, energy)')
    con.commit()

    n = con.execute('SELECT COUNT(*) FROM coincidences').fetchone()[0]
    co60 = con.execute(
        'SELECT intensity, unc_intensity FROM coincidences WHERE '
        "isotope='60CO' AND isomer='g' AND ABS(energy-1173.228)<0.1 AND "
        "ABS(coinc_energy-1332.492)<0.1 AND coinc_type='g'").fetchone()
    print(f'coincidences: {n} rows from {co_report["datasets_matched"]} '
          f'datasets')
    print(f'  60CO 1173&1332 anchor: {co60} (expect ~99.83)')
    print(f'  datasets with no matching chart state: '
          f'{len(co_report["datasets_unmatched"])} — '
          f'{co_report["datasets_unmatched"][:6]}; band-level parents '
          f'skipped: {len(co_report.get("datasets_band_parent", []))}')
    print(f'  pair energies with no gammas row to snap to: '
          f'{co_report["snap_missed"]}')
    print(f'  below-threshold directional rows dropped: '
          f'{co_report["threshold_dropped"]}; rows with no source '
          f'intensity: {co_report.get("missing_intensity", 0)}; '
          f'same-key rows aggregated: '
          f'{co_report.get("aggregated_pairs", 0)}; pairs exceeding '
          f'either single: {co_report["over_singles"]}')
    con.close()


if __name__ == '__main__':
    main()
