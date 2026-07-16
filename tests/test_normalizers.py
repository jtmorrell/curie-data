"""Tests for the shared nuclide/target normalizers.

Run with `python3 -m pytest tests/` or directly with `python3 tests/test_normalizers.py`.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts' / 'build'))
from common import curie_nuclide, curie_target


@pytest.mark.parametrize('token,expected', [
    # mass-first (IAEA reaction strings)
    ('22Na', '22NA'),
    ('48V', '48V'),
    ('94mTc', '94TCm'),
    ('102gRh', '102RHg'),
    ('99mTc', '99TCm'),
    ('117mSn', '117SNm'),
    ('52Mnm', '52MNm'),        # trailing marker, mass-first
    ('192m1Ir', '192IRm1'),
    # element-first (TALYS residual-table headers)
    ('Mn52m', '52MNm'),
    ('Mn56', '56MN'),
    ('Co57', '57CO'),
    ('Tc94m', '94TCm'),
])
def test_curie_nuclide(token, expected):
    assert curie_nuclide(token) == expected


@pytest.mark.parametrize('token', [
    'garbage', '52', 'Na', '94mTcm',  # two isomer markers
])
def test_curie_nuclide_rejects(token):
    with pytest.raises(ValueError):
        curie_nuclide(token)


@pytest.mark.parametrize('token,expected', [
    ('natTi', 'natTI'),
    ('natCu', 'natCU'),
    ('27Al', '27AL'),
    ('60Ni', '60NI'),
    ('226Ra', '226RA'),
])
def test_curie_target(token, expected):
    assert curie_target(token) == expected


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-q']))
