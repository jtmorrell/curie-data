"""Upload a data generation to its GitHub release tags.

Creates one DRAFT release per tag (drafts publish no git tag and are
invisible to the fetch URLs until published) and uploads that tag's assets:
the whole files, SHA256SUMS and data_manifest.toml on the generation tag,
and each shard directory on its own <generation>-<library> tag. The run is
resumable: an asset already present with the right size is skipped, one
with the wrong size is replaced, so an interrupted or rate-limited run is
finished by running again. Verification lists every release's assets and
compares names against the curie registry; publishing (which creates the
tags and makes the download URLs live) is a separate explicit step.

Everything uses the REST API only (the maintainer token cannot use
GraphQL); the token comes from GITHUB_TOKEN or `gh auth token`.

Usage:
    python3 release_upload.py upload [--generation v2]
    python3 release_upload.py verify [--generation v2]
    python3 release_upload.py publish [--generation v2]
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


def _token():
    tok = os.environ.get('GITHUB_TOKEN')
    if not tok:
        tok = subprocess.run(['gh', 'auth', 'token'], capture_output=True,
                             text=True, check=True).stdout.strip()
    return tok


def _request(url, method='GET', data=None, content_type='application/json',
             retries=4):
    for attempt in range(retries + 1):
        req = urllib.request.Request(url, method=method, data=data, headers={
            'Authorization': 'Bearer ' + _TOKEN,
            'Accept': 'application/vnd.github+json',
            'X-GitHub-Api-Version': '2022-11-28',
            'Content-Type': content_type})
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                body = resp.read()
                return json.loads(body) if body else None
        except urllib.error.HTTPError as e:
            # 422 on an asset upload means the name already exists (a
            # concurrent or interrupted upload won the race): caller skips
            if e.code == 422 and '/assets?' in url:
                return None
            if e.code in (500, 502, 503, 504) and attempt < retries:
                time.sleep(2 ** attempt)
                continue
            sys.exit('{} {} -> HTTP {}: {}'.format(
                method, url, e.code, e.read().decode(errors='replace')[:300]))
        except (urllib.error.URLError, TimeoutError):
            if attempt < retries:
                time.sleep(2 ** attempt)
                continue
            raise


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


def _plan(generation):
    """[(tag, [paths])] for the generation: whole files + per-library shards."""
    build = ROOT / 'build'
    reg = json.loads(REGISTRY.read_text())
    whole = [build / f for f in sorted(reg['files'])]
    missing = [str(p) for p in whole if not p.is_file()]
    if missing:
        sys.exit('whole files missing from build/: {}'.format(missing))
    sums = build / 'SHA256SUMS'
    sums.write_text(''.join(
        '{}  {}\n'.format(hashlib.sha256(p.read_bytes()).hexdigest(), p.name)
        for p in whole))
    plan = [(generation, whole + [sums, ROOT / 'data_manifest.toml'])]
    for lib, group in sorted(reg['shards'].items()):
        tag = group['base_url'].rstrip('/').rsplit('/', 1)[-1]
        shard_dir = ROOT / 'shards-v2' / lib
        paths = [shard_dir / f for f in sorted(group['files'])]
        missing = [str(p) for p in paths if not p.is_file()]
        if missing:
            sys.exit('{}: shards missing on disk: {}...'.format(tag, missing[:3]))
        plan.append((tag, paths))
    return plan


def upload(generation):
    releases = _releases()
    for tag, paths in _plan(generation):
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
        n_up, n_skip = 0, 0
        for p in paths:
            a = have.get(p.name)
            if a is not None:
                if a['size'] == p.stat().st_size:
                    n_skip += 1
                    continue
                _request(f'{API}/repos/{REPO}/releases/assets/{a["id"]}',
                         method='DELETE')
            _request(f'{UPLOADS}/repos/{REPO}/releases/{rel["id"]}/assets'
                     f'?name={p.name}', method='POST', data=p.read_bytes(),
                     content_type='application/octet-stream')
            n_up += 1
            if n_up % 25 == 0:
                print('  {}: {}/{} uploaded'.format(tag, n_up, len(paths)),
                      flush=True)
        print('{}: {} uploaded, {} already present, {} total'.format(
            tag, n_up, n_skip, len(paths)), flush=True)


def verify(generation):
    releases = _releases()
    ok = True
    for tag, paths in _plan(generation):
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
        print('{} ({}): {} assets; absent {}, extra {}, size-mismatched {}'.format(
            tag, state, len(have), absent[:3] or 0, extra[:3] or 0, bad[:3] or 0))
        ok = ok and not (absent or extra or bad)
    print('VERIFY: ' + ('OK' if ok else 'INCOMPLETE'))
    return ok


def publish(generation):
    if not verify(generation):
        sys.exit('refusing to publish an incomplete upload')
    releases = _releases()
    for tag, _ in _plan(generation):
        rel = releases[tag]
        if rel['draft']:
            _request(f'{API}/repos/{REPO}/releases/{rel["id"]}', method='PATCH',
                     data=json.dumps({'draft': False}).encode())
            print('published {}'.format(tag), flush=True)
        else:
            print('{} already published'.format(tag), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('command', choices=['upload', 'verify', 'publish'])
    ap.add_argument('--generation', default='v2')
    args = ap.parse_args()
    _TOKEN = _token()
    {'upload': upload, 'verify': verify, 'publish': publish}[args.command](args.generation)
