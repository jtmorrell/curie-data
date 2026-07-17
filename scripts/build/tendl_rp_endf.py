"""Residual-production extraction from a TENDL ENDF-6 file.

Used for targets whose TALYS tables exist on no public host — the proton
niobium gap: TENDL-2025 dropped 91-95Nb from the only public proton tree,
so tendl_p_rp carries them from the TENDL-2023 ENDF-6 files (maintainer
ruling; the mixed generation is recorded in the database's _version row).

The extraction requires the s30 representation, where MT5 plus the
low-energy exclusive channels carry the complete nonelastic flux and the
MF6 MT5 residual yields sum to one heavy residual per event — TENDL's own
rp tables are exactly reproducible from it (proven by the Fe056 identity
check in the validator). The -mt representation is NON-conservative for
residuals (its exclusive channels plus MT5 do not sum to the nonelastic)
and is rejected by the balance guard.

Per product: sigma(E) = sum over exclusive MTs whose emitted-particle
bookkeeping lands on the product of the MF3 cross section (per-level MF10
sections where the residual has metastable states), plus the MF6 MT5
yield times the MT5 cross section per LIP. State identifiers differ by
channel — MF10's LFS is a level INDEX (93Mo's isomer is LFS 16, mirroring
the TALYS .L## files) while MF6's LIP is an isomer ordinal — so each
channel's states are ranked into the shared g/m1/m2 ordinals; a channel
that lists excited levels without the ground state would break the
ranking and stops the build.

Fidelity limit of the source (documented by the Fe056 identity check):
in the few MeV below the exclusive-to-MT5 handoff (~24-30 MeV) the ENDF-6
file carries slightly less than the TALYS tables — production through
channels the exclusive MT set cannot encode (such as (p,3p)) only enters
once MT5 takes over. The observed shortfall is bounded by ~2% of the
affected product's peak (worst cases: 1.4% relative on 56Co, and 54Cr's
below-handoff tail at 1.7% of its peak); the identity check's combined
relative-and-absolute tolerance reflects exactly this.
"""
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from common import ARCHIVE_ROOT, ELEMENTS

# emitted (Z, A) per exclusive MT, incident-proton sublibrary
EMIT = {4: (0, 1), 11: (1, 4), 16: (0, 2), 17: (0, 3), 22: (2, 5),
        23: (6, 13), 24: (2, 6), 25: (2, 7), 28: (1, 2), 29: (4, 9),
        30: (4, 10), 32: (1, 3), 33: (1, 4), 34: (2, 4), 35: (5, 11),
        36: (5, 12), 37: (0, 4), 41: (1, 3), 42: (1, 4), 44: (2, 3),
        45: (3, 6), 102: (0, 0), 103: (1, 1), 104: (1, 2), 105: (1, 3),
        106: (2, 3), 107: (2, 4), 108: (4, 8), 109: (6, 12), 111: (2, 2),
        112: (3, 5), 113: (5, 11), 114: (5, 10), 115: (2, 3),
        116: (2, 4), 117: (3, 6)}
NON_CHANNEL_MTS = {1, 2, 3, 5}  # totals/elastic/sums; MT5 handled via MF6

# Evaluation backbone [MeV]: the source grids are coarse (a few tens of
# points) and encode the ~30 MeV exclusive->MT5 handoff as duplicate-
# abscissa steps that a deduplicated union grid would smear, so every
# contribution is evaluated (through its own interpolation laws) on this
# fine grid unioned with the source points. Quarter-MeV steps resolve the
# handoff four times finer than the TALYS tables' own 1 MeV grid.
BACKBONE_MEV = np.arange(0.25, 200.001, 0.25)


def _eval(tab, e_ev):
    """Evaluate an endf Tabulated1D at each energy (its own interpolation
    laws apply within each contribution's source grid)."""
    return np.array([tab(e) for e in e_ev], float)


def endf_targets(mirror_src):
    """(curie_target, table_name, path) for each archived ENDF-6 target
    not already covered by the table mirror under `mirror_src`."""
    root = ARCHIVE_ROOT / 'tendl-2023' / 'proton_endf'
    out = []
    for p in sorted(root.rglob('p-*.tendl')):
        m = re.fullmatch(r'p-([A-Z][a-z]?)(\d{3})\.tendl', p.name)
        el, a = m.group(1), int(m.group(2))
        if (mirror_src / el / f'{el}{a:03d}' / 'tables' / 'residual').is_dir():
            print(f'  {p.name}: table mirror now covers {el}{a} — ENDF-6 '
                  f'ingest skipped')
            continue
        out.append((f'{a}{el.upper()}', f'{el.upper()}_{a}', p))
    return out


def extract_products(path, e_max_mev=200.0):
    """{product_name: rows[(E [MeV], sigma [mb])]} from an s30-style file."""
    import endf
    sd = endf.Material(str(path)).section_data

    s3 = sd[(3, 3)]['sigma']
    s5 = sd[(3, 5)]['sigma']
    z_t, a_t = divmod(int(round(sd[(3, 5)]['ZA'])), 1000)
    zc, ac = z_t + 1, a_t + 1

    # s30-representation guard: at the top of the grid the exclusive
    # channels have handed everything to MT5, so excl+MT5 == nonelastic
    e_top = min(s5.x[-1], e_max_mev * 1e6)
    excl_top = sum(float(sd[(3, mt)]['sigma'](e_top))
                   for mf, mt in sd if mf == 3 and mt not in NON_CHANNEL_MTS)
    balance = (excl_top + float(s5(e_top))) / float(s3(e_top))
    if not 0.99 < balance < 1.01:
        raise SystemExit(
            f'{path}: exclusives+MT5 = {balance:.3f} x nonelastic at '
            f'{e_top / 1e6:.0f} MeV — not an s30-representation file '
            f'(the -mt variant cannot reconstruct residual production)')

    # (z, a) -> list of (state_ordinal | None, source grid [eV], mb_fn(E_eV))
    contribs = {}

    for p in sd[(6, 5)]['products']:
        z, a = divmod(int(p['ZAP']), 1000)
        if z <= 2 or a == 0:
            continue                      # ejectiles and photons, not residuals
        y = p['y_i']
        grid = np.concatenate([np.asarray(y.x, float), np.asarray(s5.x, float)])
        contribs.setdefault((z, a), []).append(
            (int(p['LIP']), grid,
             lambda e, y=y: _eval(s5, e) * _eval(y, e) * 1e3))

    for (mf, mt), sec in sorted(sd.items()):
        if mf != 3 or mt in NON_CHANNEL_MTS:
            continue
        if mt not in EMIT:
            raise SystemExit(f'{path}: unmapped exclusive channel MT{mt}')
        dz, da = EMIT[mt]
        z, a = zc - dz, ac - da
        if not 2 < z < len(ELEMENTS):
            continue
        if (10, mt) in sd:
            levels = sd[(10, mt)]['levels']
            lfs = [int(lv['LFS']) for lv in levels]
            if sorted(lfs)[0] != 0:
                raise SystemExit(f'{path}: MF10 MT{mt} lists excited levels '
                                 f'without the ground state — cannot rank')
            for rank, lv in enumerate(sorted(levels,
                                             key=lambda l: int(l['LFS']))):
                sig = lv['sigma']
                contribs.setdefault((z, a), []).append(
                    (rank, np.asarray(sig.x, float),
                     lambda e, sig=sig: _eval(sig, e) * 1e3))
        else:
            # no MF8/MF10: the residual has no metastable state and the
            # whole channel is its ground-state production
            sig = sec['sigma']
            contribs.setdefault((z, a), []).append(
                (None, np.asarray(sig.x, float),
                 lambda e, sig=sig: _eval(sig, e) * 1e3))

    products = {}
    for (z, a), parts in contribs.items():
        split = any(o not in (None, 0) for o, *_ in parts)
        by_state = {}
        for o, grid, fn in parts:
            if o is None:
                if split and fn(grid).max() > 1e-6:
                    raise SystemExit(
                        f'{path}: unsplit channel feeds isomeric product '
                        f'Z={z} A={a} — state attribution impossible')
                o = 0
            by_state.setdefault(o, []).append((grid, fn))
        for o, series in by_state.items():
            e = np.unique(np.concatenate(
                [BACKBONE_MEV * 1e6] + [g for g, _ in series]))
            e = e[e <= e_max_mev * 1e6]
            tot = np.zeros_like(e)
            for _, fn in series:
                tot += fn(e)
            name = f'{a}{ELEMENTS[z].upper()}' + ('g' if o == 0 else f'm{o}')
            products[name] = list(zip((e / 1e6).tolist(), tot.tolist()))
    return products
