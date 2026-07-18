"""Upload a data generation to its GitHub release tags.

Creates one DRAFT release per tag (drafts publish no git tag and are
invisible to the fetch URLs until published) and uploads that tag's assets:
the whole files, SHA256SUMS and data_manifest.toml on the generation tag,
and each shard directory on its own <generation>-<library> tag. The
generation and every tag name derive from the curie registry (base_url and
the per-group base_urls), so the uploaded layout can never disagree with
what curie will fetch. The run is resumable: an asset already present with
the right size is skipped, one with the wrong size is replaced, and the
two metadata files (SHA256SUMS, data_manifest.toml) are always replaced —
their size barely moves when their content changes, so size is no proxy
for them. API rate-limit exhaustion sleeps until the advertised reset.

`verify` first proves the LOCAL files match the registry checksums
byte-for-byte (whole files and every shard), then that every release
holds exactly the planned asset names at the local sizes; since uploads
are single atomic POSTs, hash-verified-local + name/size-verified-remote
pins the published set. `publish` (which creates the tags and makes the
download URLs live) refuses to run on anything but a clean verify.

Everything uses the REST API only (the maintainer token cannot use
GraphQL); the token comes from GITHUB_TOKEN or `gh auth token`.

Usage:
    python3 release_upload.py upload
    python3 release_upload.py verify
    python3 release_upload.py publish
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

REPO = 'jtmorrell/curie-data'
API = 'https://api.github.com'
UPLOADS = 'https://uploads.github.com'
ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT.parent / 'curie' / 'curie' / 'data_registry.json'
METADATA = ('SHA256SUMS', 'data_manifest.toml')  # size-stable: always replace
# per-asset SHA256s recorded at upload time: verify compares them against the
# current local files, catching an in-place content change that kept the byte
# count (twice observed for decay.db) which the remote size check cannot see
LEDGER = ROOT / 'build' / 'uploaded_hashes.json'

_TOKEN = None


def _token():
    global _TOKEN
    if _TOKEN is None:
        _TOKEN = os.environ.get('GITHUB_TOKEN') or subprocess.run(
            ['gh', 'auth', 'token'], capture_output=True, text=True,
            check=True).stdout.strip()
    return _TOKEN


def _request(url, method='GET', data=None, content_type='application/json',
             retries=4):
    attempt, rate_waits = 0, 0
    while True:
        req = urllib.request.Request(url, method=method, data=data, headers={
            'Authorization': 'Bearer ' + _token(),
            'Accept': 'application/vnd.github+json',
            'X-GitHub-Api-Version': '2022-11-28',
            'Content-Type': content_type})
        try:
            with urllib.request.urlopen(req, timeout=600) as resp:
                body = resp.read()
                return json.loads(body) if body else None
        except urllib.error.HTTPError as e:
            # 422 on an asset upload means the name already exists (a
            # concurrent or interrupted upload won the race): caller skips.
            # 404 on an asset DELETE means it is already gone.
            if e.code == 422 and '/assets?' in url:
                return None
            if e.code == 404 and method == 'DELETE':
                return None
            # a full generation is ~2,900 requests against the 5,000/hour
            # limit: on exhaustion, sleep until the advertised reset and
            # keep going (rate waits do not consume error retries)
            retry_after = e.headers.get('Retry-After')
            if (e.code in (403, 429)
                    and (e.headers.get('X-RateLimit-Remaining') == '0' or retry_after)
                    and rate_waits < 4):
                rate_waits += 1
                # secondary (abuse) limits advertise Retry-After instead of
                # a reset timestamp
                wait = (min(int(retry_after) + 5, 4000) if retry_after else
                        min(max(int(e.headers.get('X-RateLimit-Reset', 0)) -
                                time.time() + 30, 60), 4000))
                print('rate limit reached; resuming in {:.0f} min'.format(wait / 60),
                      flush=True)
                time.sleep(wait)
                continue
            if e.code in (500, 502, 503, 504) and attempt < retries:
                attempt += 1
                time.sleep(2 ** attempt)
                continue
            sys.exit('{} {} -> HTTP {}: {}'.format(
                method, url, e.code, e.read().decode(errors='replace')[:300]))
        except (urllib.error.URLError, TimeoutError):
            if attempt >= retries:
                raise
            attempt += 1
            time.sleep(2 ** attempt)


def _releases():
    """{tag_name: release dict} including drafts (which /releases/tags/
    cannot see)."""
    out, page = {}, 1
    while True:
        batch = _request(f'{API}/repos/{REPO}/releases?per_page=100&page={page}')
        if not batch:
            return out
        for r in batch:
            out[r['tag_name']] = r
        page += 1


def _assets(release_id):
    out, page = {}, 1
    while True:
        batch = _request(f'{API}/repos/{REPO}/releases/{release_id}/assets'
                         f'?per_page=100&page={page}')
        if not batch:
            return out
        for a in batch:
            out[a['name']] = a
        page += 1


def _plan(write_sums=False):
    """(generation, [(tag, [paths])]) from the registry: whole files +
    metadata on the generation tag, each shard group on its own tag. Only
    an upload regenerates SHA256SUMS; verify/publish read what upload
    wrote."""
    build = ROOT / 'build'
    reg = json.loads(REGISTRY.read_text())
    generation = reg['base_url'].rstrip('/').rsplit('/', 1)[-1]
    whole = [build / f for f in sorted(reg['files'])]
    missing = [str(p) for p in whole if not p.is_file()]
    if missing:
        sys.exit('whole files missing from build/: {}'.format(missing))
    sums = build / 'SHA256SUMS'
    if write_sums:
        sums.write_text(''.join(
            '{}  {}\n'.format(hashlib.sha256(p.read_bytes()).hexdigest(), p.name)
            for p in whole))
    manifest = ROOT / 'data_manifest.toml'
    for p in (sums, manifest):
        if not p.is_file():
            sys.exit('{} missing (run upload first)'.format(p))
    plan = [(generation, whole + [sums, manifest])]
    for lib, group in sorted(reg['shards'].items()):
        tag = group['base_url'].rstrip('/').rsplit('/', 1)[-1]
        shard_dir = ROOT / 'shards-{}'.format(generation) / lib
        paths = [shard_dir / f for f in sorted(group['files'])]
        missing = [str(p) for p in paths if not p.is_file()]
        if missing:
            sys.exit('{}: shards missing on disk: {}...'.format(tag, missing[:3]))
        plan.append((tag, paths))
    return generation, plan


def upload():
    generation, plan = _plan(write_sums=True)
    releases = _releases()
    for tag, paths in plan:
        rel = releases.get(tag)
        if rel is None:
            rel = _request(f'{API}/repos/{REPO}/releases', method='POST',
                           data=json.dumps({
                               'tag_name': tag, 'name': tag, 'draft': True,
                               'body': 'curie nuclear-data release, generation '
                                       '{}. Fetched by curie via '
                                       'curie/data_registry.json.'.format(generation),
                           }).encode())
            print('created DRAFT release {}'.format(tag), flush=True)
        have = _assets(rel['id'])
        ledger = json.loads(LEDGER.read_text()) if LEDGER.exists() else {}
        rec = ledger.setdefault(tag, {})
        n_up, n_skip = 0, 0
        for p in paths:
            data = p.read_bytes()
            sha = hashlib.sha256(data).hexdigest()
            a = have.get(p.name)
            if a is not None:
                if (a['size'] == len(data) and p.name not in METADATA
                        and rec.get(p.name) in (None, sha)):
                    rec[p.name] = rec.get(p.name, sha)
                    n_skip += 1
                    continue
                _request(f'{API}/repos/{REPO}/releases/assets/{a["id"]}',
                         method='DELETE')
            _request(f'{UPLOADS}/repos/{REPO}/releases/{rel["id"]}/assets'
                     f'?name={p.name}', method='POST', data=data,
                     content_type='application/octet-stream')
            rec[p.name] = sha
            n_up += 1
            if n_up % 25 == 0:
                print('  {}: {}/{} uploaded'.format(tag, n_up, len(paths)),
                      flush=True)
        print('{}: {} uploaded, {} already present, {} total'.format(
            tag, n_up, n_skip, len(paths)), flush=True)
        LEDGER.write_text(json.dumps(ledger, indent=1, sort_keys=True))


def verify():
    generation, plan = _plan()
    reg = json.loads(REGISTRY.read_text())
    # local disk must equal the registry checksums byte-for-byte: whole
    # files and every shard (the uploads are atomic single POSTs, so
    # hash-verified local + name/size-verified remote pins the release)
    bad_hash = []
    for f, sha in reg['files'].items():
        p = ROOT / 'build' / f
        if hashlib.sha256(p.read_bytes()).hexdigest() != sha:
            bad_hash.append(f)
    for lib, group in reg['shards'].items():
        for f, sha in group['files'].items():
            p = ROOT / 'shards-{}'.format(generation) / lib / f
            if hashlib.sha256(p.read_bytes()).hexdigest() != sha:
                bad_hash.append(f)
    n_checked = len(reg['files']) + sum(len(g['files']) for g in reg['shards'].values())
    print('local files vs registry checksums: {} checked, {} mismatched{}'.format(
        n_checked, len(bad_hash), ' ' + str(bad_hash[:5]) if bad_hash else ''))
    ok = not bad_hash
    # an asset whose local content changed since its upload keeps its byte
    # count more often than not; the remote size check cannot see that
    if LEDGER.exists():
        ledger = json.loads(LEDGER.read_text())
        stale = []
        for tag, paths in plan:
            rec = ledger.get(tag, {})
            for pth in paths:
                if pth.name in rec and rec[pth.name] != hashlib.sha256(
                        pth.read_bytes()).hexdigest():
                    stale.append('{}/{}'.format(tag, pth.name))
        print('local files changed since their upload: {}{}'.format(
            len(stale), ' ' + str(stale[:5]) if stale else ''))
        ok = ok and not stale

    releases = _releases()
    for tag, paths in plan:
        rel = releases.get(tag)
        if rel is None:
            print('{}: MISSING release'.format(tag))
            ok = False
            continue
        have = _assets(rel['id'])
        wanted = {p.name: p.stat().st_size for p in paths}
        absent = sorted(set(wanted) - set(have))
        extra = sorted(set(have) - set(wanted))
        bad = sorted(n for n in wanted if n in have
                     and have[n]['size'] != wanted[n])
        state = 'draft' if rel['draft'] else 'published'
        def _fmt(lst):
            return '{} {}'.format(len(lst), lst[:3]) if lst else '0'
        print('{} ({}): {} assets; absent {}, extra {}, size-mismatched {}'.format(
            tag, state, len(have), _fmt(absent), _fmt(extra), _fmt(bad)))
        ok = ok and not (absent or extra or bad)
    print('VERIFY: ' + ('OK' if ok else
          'INCOMPLETE — fix the mismatches above, re-run upload, then verify'))
    return ok


def publish():
    if not verify():
        sys.exit('refusing to publish an incomplete or mismatched upload')
    _, plan = _plan()
    releases = _releases()
    for tag, _ in plan:
        rel = releases.get(tag)
        if rel is None:
            sys.exit('release {} vanished between verify and publish'.format(tag))
        if rel['draft']:
            _request(f'{API}/repos/{REPO}/releases/{rel["id"]}', method='PATCH',
                     data=json.dumps({'draft': False}).encode())
            print('published {}'.format(tag), flush=True)
        else:
            print('{} already published'.format(tag), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('command', choices=['upload', 'verify', 'publish'])
    args = ap.parse_args()
    {'upload': upload, 'verify': verify, 'publish': publish}[args.command]()
