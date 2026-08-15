"""Tests for the half-life-aware coincidence pair cut in decay_coinc.

A gamma separated from its partner photon by a long-lived intermediate level
is not truly coincident: the cut drops any pair whose cascade transits a level
whose half-life exceeds ISOMER_HALFLIFE_CUT_S. Fixtures mirror the paceENSDF
JSON structure (decaySingles level records, GG indexed cascade sequences, GX
single-hop gamma legs) rather than the shipped data.

Run with `python3 -m pytest tests/` or directly with
`python3 tests/test_decay_coinc.py`.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts' / 'build'))
from decay_coinc import (ISOMER_HALFLIFE_CUT_S, build_singles_index,
                         cascade_transitions, transits_delaying_isomer)


def _single(i_init, i_final, isomer=False, hl=None):
    """A decaySingles level record with the fields the cut reads."""
    return {'levelIndexInitial': i_init, 'levelIndexFinal': i_final,
            'levelIsomer': isomer, 'levelHalfLifeConverted': hl}


# A daughter level scheme modeled on 113In following 113Sn EC/beta+ decay:
# the 391.698 keV level (index 1 -> 0) is the 5968.56 s 113mIn isomer; the
# 1029.73 keV level (index 3) and 646.83 keV level (index 2) are prompt.
SINGLES = [
    _single(1, 0, isomer=True, hl=5968.56),   # 391.698 keV isomer
    _single(2, 0),                             # prompt 646.83 keV
    _single(2, 1),                             # prompt 255.134 keV
    _single(3, 2, isomer=True, hl=3.3e-10),    # sub-ns isomer (below the cut)
    _single(3, 1),                             # prompt 638.03 keV
    _single(4, 0, isomer=True, hl=None),       # isomer, unknown half-life
]
INDEX = build_singles_index(SINGLES)


@pytest.mark.parametrize('seq,expect_cut,expect_unknown', [
    # GG mid-cascade isomer: 255.134/391.698 cascade [2, 1, 0] transits the
    # 391.698 keV isomer on its (1 -> 0) leg -> cut. The gate level (2) is
    # prompt, so the gate's own half-life field would miss this.
    ([2, 1, 0], True, False),
    # GG prompt cascade: 646.83/255.134 [2, 1] then... use [3, 1] (638.03)
    # straight to the 391.698 level's floor without leaving an isomer.
    ([3, 1], False, False),
    # GG cascade transiting only a sub-nanosecond isomer (3 -> 2) stays.
    ([3, 2], False, False),
    # GG cascade transiting the unknown-half-life isomer (4 -> 0) -> cut,
    # flagged unknown.
    ([4, 0], True, True),
])
def test_gg_cascade_cut(seq, expect_cut, expect_unknown):
    transitions = cascade_transitions(seq)
    assert transits_delaying_isomer(transitions, INDEX) == (expect_cut,
                                                            expect_unknown)


@pytest.mark.parametrize('gamma_leg,expect_cut,expect_unknown', [
    # GX blocked: the 391.698 keV gamma leg (1 -> 0) leaves the isomer.
    ((1, 0), True, False),
    # GX passing: the 255.134 keV gamma leg (2 -> 1) is prompt.
    ((2, 1), False, False),
    # GX blocked by an unknown-half-life isomer level.
    ((4, 0), True, True),
])
def test_gx_single_hop_cut(gamma_leg, expect_cut, expect_unknown):
    assert transits_delaying_isomer([gamma_leg], INDEX) == (expect_cut,
                                                            expect_unknown)


def test_threshold_boundary():
    """A level exactly at the cut is kept; one above it is dropped."""
    idx = build_singles_index([
        _single(1, 0, isomer=True, hl=ISOMER_HALFLIFE_CUT_S),
        _single(2, 0, isomer=True, hl=ISOMER_HALFLIFE_CUT_S * 10),
    ])
    assert transits_delaying_isomer([(1, 0)], idx) == (False, False)
    assert transits_delaying_isomer([(2, 0)], idx) == (True, False)


def test_non_isomer_never_cuts():
    """A long half-life on a non-isomer level never cuts."""
    idx = build_singles_index([_single(1, 0, isomer=False, hl=1e9)])
    assert transits_delaying_isomer([(1, 0)], idx) == (False, False)


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-q']))
