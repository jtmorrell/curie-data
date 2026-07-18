"""Mirror the TENDL-2025 tabulated cross sections for the curated target set.

Per beam (neutron, proton, deuteron, alpha) and target: the residual-
production tables (`tables/residual/rp*.tot` and isomer `.L##` files) and,
for the neutron beam only, the exclusive-channel tables (`tables/xs/*.tot`
plus the total/elastic/nonelastic/fission/*prod specials). File URLs are
enumerated from the per-target HTML index pages, whose references are
plain text rather than anchors.

This mirror is ~half a million small files served by a static university
archive, so it runs with a bounded worker pool and connection keep-alive
instead of the one-request-per-second throttle used for the IAEA's
dynamic services. Provenance is recorded as one digest per (beam, target)
in the shared manifest — the combined SHA256 over the target's sorted
per-file digests, which are themselves written to a MANIFEST.sha256 in
the target's archive directory — keeping the committed manifest at one
entry per target rather than one per file. On resume, existing files are
re-hashed against that MANIFEST.sha256 and re-fetched on mismatch, and
the shared manifest gains a new entry only when a target's digest
actually changed.

Usage:
    python3 tendl.py <targets-file> [beam ...]     # beams default to all
"""
import hashlib
import re
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from core import ARCHIVE_ROOT, _record, _recorded_sha

BASE = 'https://tendl.imperial.ac.uk/tendl_2025'
# The Imperial tree is the only public host: the PSI host is gone, and the
# NEA Data Bank copy sits behind a GitLab login whose file paths redirect
# (3xx) to a sign-in page — http_get raises on any non-200 status, so
# nothing from it can enter the archive. Imperial's known gaps (its proton
# sublibrary lacks the Nb isotopes) are reported as not-provided.
UA = 'curie-data-pipeline/2.0 (github.com/jtmorrell/curie-data)'
WORKERS = 12
BEAMS = ['neutron', 'proton', 'deuteron', 'alpha']

ELEMENTS = ['n','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',
            'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',
            'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr',
            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I',
            'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
            'Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au',
            'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U',
            'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db',
            'Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts']
_BY_UPPER = {e.upper(): e for e in ELEMENTS}


_local = threading.local()


def http_get(url, timeout=60, _retry=True):
    """GET over a per-thread persistent connection (the mirror is hundreds
    of thousands of small files; per-request TLS setup would dominate)."""
    import http.client
    from urllib.parse import urlsplit
    parts = urlsplit(url)
    conn = getattr(_local, 'conn', None)
    if conn is None or getattr(_local, 'host', None) != parts.netloc:
        if conn is not None:
            conn.close()
        conn = http.client.HTTPSConnection(parts.netloc, timeout=timeout)
        _local.conn, _local.host = conn, parts.netloc
    try:
        conn.request('GET', parts.path, headers={'User-Agent': UA})
        resp = conn.getresponse()
        data = resp.read()
    except (http.client.HTTPException, OSError):
        _local.conn = None
        if _retry:
            return http_get(url, timeout, _retry=False)  # one fresh-connection retry
        raise
    if resp.status == 404:
        raise FileNotFoundError(url)
    if resp.status != 200:
        raise OSError(f'HTTP {resp.status} for {url}')
    expected = resp.headers.get('Content-Length')
    if expected is not None and len(data) != int(expected):
        raise OSError(f'short body for {url}')
    return data


def tendl_paths(curie_target):
    """'152EU' / '110AGm' / '178HFm2' -> (element, A, TENDL isomer letter,
    padded nuclide token): ('Eu', 152, '', 'Eu152') / ('Ag', 110, 'm',
    'Ag110m') / ('Hf', 178, 'n', 'Hf178n'). TALYS isomer-target letters
    run m, n, o, ... in step with the curie ordinals m(=m1), m2, m3."""
    m = re.fullmatch(r'(\d+)([A-Z]{1,2})(m(\d+)?)?', curie_target)
    a, el, iso, k = m.groups()
    letter = chr(ord('m') + int(k or 1) - 1) if iso else ''
    element = _BY_UPPER[el]
    return element, int(a), letter, f'{element}{int(a):03d}{letter}'


def index_pages(beam, curie_target):
    """The target's HTML index-page URLs. Page URLs carry the UNPADDED mass
    number while file paths use the zero-padded nuclide token — both
    conventions are TENDL's own (pinned by tests/test_tendl_fetch.py)."""
    element, a, iso, _ = tendl_paths(curie_target)
    cap = beam.capitalize()
    pages = [f'{BASE}/{beam}_html/{element}/{cap}{element}{a}{iso}residual.html']
    if beam == 'neutron':
        pages.append(pages[0].replace('residual.html', 'xs.html'))
    return pages


def file_urls(beam, curie_target):
    """Enumerate this target's table files from the HTML index pages.
    Returns (urls, missing_flag). A 404'd index page sets the missing flag
    (the target, or one of its pages, is genuinely absent upstream); any
    other failure raises, so a transient enumeration error surfaces in the
    failed tally instead of masquerading as 'not provided'."""
    urls, missing = set(), False
    for page in index_pages(beam, curie_target):
        try:
            html = http_get(page).decode('utf-8', 'replace')
        except FileNotFoundError:
            missing = True
            continue
        urls |= {'https:' + u if u.startswith('//') else u for u in
                 re.findall(r'[\w/.\-:]+/tables/(?:residual|xs)/[\w.]+', html)}
    return sorted(urls), missing


def fetch_target(beam, curie_target):
    """Mirror one (beam, target); returns (target, n_files, digest|None, note)."""
    urls, missing = file_urls(beam, curie_target)
    if not urls:
        return curie_target, 0, None, 'not provided'
    element, _, _, nuc = tendl_paths(curie_target)
    dest_dir = ARCHIVE_ROOT / 'tendl-2025' / f'{beam}_files' / element / nuc / 'tables'
    listing_file = dest_dir / 'MANIFEST.sha256'
    recorded = {}
    if listing_file.exists():
        for line in listing_file.read_text().splitlines():
            sha, _, rel = line.partition('  ')
            recorded[rel] = sha
    shas = {}
    for url in urls:
        rel = url.split('/tables/')[-1]
        dest = dest_dir / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and dest.stat().st_size > 0:
            sha = hashlib.sha256(dest.read_bytes()).hexdigest()
            if recorded.get(rel, sha) != sha:
                dest.unlink()  # corrupt since recorded: re-fetch below
        if not (dest.exists() and dest.stat().st_size > 0):
            data = http_get(url)
            part = dest.with_name(dest.name + '.part')
            part.write_bytes(data)
            part.replace(dest)
            sha = hashlib.sha256(data).hexdigest()
        shas[rel] = sha
    listing = ''.join(f'{shas[k]}  {k}\n' for k in sorted(shas))
    listing_file.write_text(listing)
    digest = hashlib.sha256(listing.encode()).hexdigest()
    source = f'tendl-2025/{beam}'
    if _recorded_sha(source, nuc) != digest:  # unchanged targets append nothing
        _record(source, nuc, f'{BASE}/{beam}_files/{element}/{nuc}/tables/',
                digest, n_files=len(shas))
    return curie_target, len(shas), digest, 'partial (some index pages missing)' if missing else 'ok'


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__.strip())
    targets = [t.strip() for t in open(sys.argv[1]) if t.strip()]
    beams = sys.argv[2:] or BEAMS
    unknown = [b for b in beams if b not in BEAMS]
    if unknown:
        sys.exit(f'unknown beam(s): {" ".join(unknown)} — choose from: {" ".join(BEAMS)}')
    for beam in beams:
        print(f'=== {beam}: {len(targets)} targets ===', flush=True)
        done, processed, not_provided, partial, failed = 0, 0, [], [], []
        with ThreadPoolExecutor(WORKERS) as pool:
            futures = {pool.submit(fetch_target, beam, t): t for t in targets}
            for fut in as_completed(futures):
                t = futures[fut]
                processed += 1
                try:
                    tgt, n, digest, note = fut.result()
                except Exception as e:
                    failed.append((t, repr(e)))
                    print(f'  {t}: FAILED {e}', flush=True)
                    continue
                if note == 'not provided':
                    not_provided.append(tgt)
                    continue
                if note.startswith('partial'):
                    partial.append(tgt)
                done += 1
                if processed % 25 == 0:
                    print(f'  ... {processed}/{len(targets)}', flush=True)
        print(f'{beam}: {done} targets mirrored ({len(not_provided)} not provided '
              f'by TENDL, {len(partial)} partial, {len(failed)} failed)', flush=True)
        if not_provided:
            print('  not provided:', ' '.join(sorted(not_provided)), flush=True)
        for f_ in failed:
            print('  FAILED:', *f_, flush=True)


if __name__ == '__main__':
    main()
