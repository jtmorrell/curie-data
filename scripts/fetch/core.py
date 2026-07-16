"""Shared fetch layer for the data-source mirror.

All network access in this repository goes through `fetch()`: throttled,
with a per-host User-Agent policy and a SHA256 provenance record appended
to `sources_manifest.toml` at the repository root. Downloads are atomic
(written to a .part file, recorded, then renamed into place), verified
against the response's Content-Length, and the resume check re-verifies an
existing file against its recorded SHA256 — a file that fails verification
is re-fetched. The manifest is append-only; when a file is re-fetched the
newest entry for a (source, name) pair is the current one.

Builders never touch the network: they read exclusively from the local
source archive (default `~/Documents/nucleardata/data_sources/v2/`,
override with CURIE_SOURCE_ARCHIVE), so a build is reproducible offline
from an existing archive.
"""
import hashlib
import os
import time
import tomllib
import urllib.request
from datetime import date
from pathlib import Path

ARCHIVE_ROOT = Path(os.environ.get(
    'CURIE_SOURCE_ARCHIVE',
    Path.home() / 'Documents/nucleardata/data_sources/v2'))

REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST = REPO_ROOT / 'sources_manifest.toml'

# The IAEA services intermittently reject clients without a browser-like
# User-Agent; the LiveChart API guide's own workaround is 'Livechart/1.0'.
_UA_DEFAULT = 'curie-data-pipeline/2.0 (github.com/jtmorrell/curie-data)'
_UA_BY_HOST = {
    'nds.iaea.org': 'Livechart/1.0',
    'www-nds.iaea.org': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36',
}

THROTTLE_S = 1.0
_last_request = [0.0]


def _throttle():
    wait = _last_request[0] + THROTTLE_S - time.monotonic()
    if wait > 0:
        time.sleep(wait)
    _last_request[0] = time.monotonic()


_recorded = None  # (source, name) -> sha256, newest entry wins; lazy-loaded


def _recorded_sha(source, name):
    global _recorded
    if _recorded is None:
        _recorded = {}
        if MANIFEST.exists():
            with open(MANIFEST, 'rb') as f:
                for entry in tomllib.load(f).get('file', []):
                    _recorded[(entry['source'], entry['name'])] = entry['sha256']
    return _recorded.get((source, name))


def _record(source, name, url, sha, size):
    entry = (
        f'\n[[file]]\n'
        f'source = "{source}"\n'
        f'name = "{name}"\n'
        f'url = "{url}"\n'
        f'sha256 = "{sha}"\n'
        f'size = {size}\n'
        f'fetched = "{date.today().isoformat()}"\n'
    )
    with open(MANIFEST, 'a') as f:
        f.write(entry)
    if _recorded is not None:
        _recorded[(source, name)] = sha


def fetch(url, source, name=None, timeout=120):
    """Fetch `url` into the archive under `source`/`name`; return the path.

    An existing file is kept only if it matches its recorded SHA256 (an
    unrecorded or mismatching file is re-fetched). Raises
    urllib.error.HTTPError on failure and OSError on a short body; callers
    decide whether a missing file is fatal.
    """
    name = name or url.rsplit('/', 1)[-1]
    dest = ARCHIVE_ROOT / source / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        sha = _recorded_sha(source, name)
        if sha and hashlib.sha256(dest.read_bytes()).hexdigest() == sha:
            return dest
        dest.unlink()  # unrecorded or corrupt: re-fetch
    _throttle()
    host = url.split('/')[2]
    req = urllib.request.Request(
        url, headers={'User-Agent': _UA_BY_HOST.get(host, _UA_DEFAULT)})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = resp.read()
        expected = resp.headers.get('Content-Length')
    if expected is not None and len(data) != int(expected):
        raise OSError(f'short body for {url}: {len(data)} of {expected} bytes')
    part = dest.with_name(dest.name + '.part')
    part.write_bytes(data)
    _record(source, name, url, hashlib.sha256(data).hexdigest(), len(data))
    os.replace(part, dest)
    return dest


def fetch_text(url, source, name=None, encoding='utf-8'):
    """fetch() and return the decoded text content."""
    return fetch(url, source, name).read_text(encoding, errors='replace')
