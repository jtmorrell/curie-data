"""Pin the TENDL mirror's URL and path conventions.

TENDL's own site uses two mass-number conventions at once: index-page URLs
carry the UNPADDED mass number (NeutronB10residual.html) while file paths
use the zero-padded nuclide token (B010). A wrong padding choice would 404
the index page and report the target as not provided, so both conventions
are pinned here without touching the network.

Run with `python3 -m pytest tests/`.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts' / 'fetch'))
from tendl import index_pages, tendl_paths


@pytest.mark.parametrize('target,expected', [
    ('10B', ('B', 10, '', 'B010')),          # single-letter element, A < 100
    ('110AGm', ('Ag', 110, 'm', 'Ag110m')),  # isomeric target
    ('183W', ('W', 183, '', 'W183')),
    ('152EU', ('Eu', 152, '', 'Eu152')),
    ('9BE', ('Be', 9, '', 'Be009')),
])
def test_tendl_paths(target, expected):
    assert tendl_paths(target) == expected


def test_index_pages_use_unpadded_mass():
    pages = index_pages('neutron', '10B')
    assert pages[0].endswith('/neutron_html/B/NeutronB10residual.html')
    assert pages[1].endswith('/neutron_html/B/NeutronB10xs.html')


def test_index_pages_isomer_and_single_page_beams():
    (page,) = index_pages('proton', '110AGm')
    assert page.endswith('/proton_html/Ag/ProtonAg110mresidual.html')
    (page,) = index_pages('alpha', '180TAm')
    assert page.endswith('/alpha_html/Ta/AlphaTa180mresidual.html')
