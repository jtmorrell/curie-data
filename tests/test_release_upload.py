"""Tests for scripts/release_upload.py and the decay stamp check of
scripts/validate.py.

Each test runs the uploader on a temporary copy of the repository layout
(build/, shards-v2/, data_manifest.toml) whose files are a few bytes each,
with the request function replaced by an in-memory release store, so no
request reaches the network and nothing under the real build/ changes. The
registry is the curie registry of the sibling checkout with its checksums
replaced by those of the small files.

Run with `python3 -m pytest tests/`.
"""
import copy
import hashlib
import json
import sqlite3
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import release_upload  # noqa: E402
import validate  # noqa: E402

CURIE_REGISTRY = ROOT.parent / 'curie' / 'curie' / 'data_registry.json'
pytestmark = pytest.mark.skipif(not CURIE_REGISTRY.is_file(),
                                reason='needs the sibling curie checkout')

API = release_upload.API
UPLOADS = release_upload.UPLOADS
REPO = release_upload.REPO


def _sha(b):
    return hashlib.sha256(b).hexdigest()


def _tag(base_url):
    return base_url.rstrip('/').rsplit('/', 1)[-1]


@pytest.fixture
def layout(tmp_path, monkeypatch):
    return make_layout(tmp_path, monkeypatch)


def make_layout(tmp_path, monkeypatch):
    """A temporary copy of the layout under the curie registry (revisions
    included), every file a few bytes and the registry's checksums those of
    the files; build/SHA256SUMS and the ledger stand for the v2 ones."""
    reg = json.loads(CURIE_REGISTRY.read_text())
    build = tmp_path / 'build'
    build.mkdir()
    for f in reg['files']:
        data = 'whole file {}\n'.format(f).encode()
        (build / f).write_bytes(data)
        reg['files'][f] = _sha(data)
    gen = _tag(reg['base_url'])
    for lib, group in reg['shards'].items():
        d = tmp_path / 'shards-{}'.format(gen) / lib
        d.mkdir(parents=True)
        for f in group['files']:
            data = 'shard {}\n'.format(f).encode()
            (d / f).write_bytes(data)
            group['files'][f] = _sha(data)
    (tmp_path / 'data_manifest.toml').write_text('[meta]\ngeneration = "v2"\n')
    sums = ''.join('{}  {}\n'.format('0' * 64, f) for f in sorted(reg['files']))
    (build / 'SHA256SUMS').write_text(sums)
    ledger = {'v2': {f: '0' * 64 for f in reg['files']},
              'v2-endf': {'endf__version.db': '0' * 64}}
    (build / 'uploaded_hashes.json').write_text(json.dumps(ledger, indent=1))
    regp = tmp_path / 'data_registry.json'
    regp.write_text(json.dumps(reg))
    monkeypatch.setattr(release_upload, 'ROOT', tmp_path)
    monkeypatch.setattr(release_upload, 'REGISTRY', regp)
    monkeypatch.setattr(release_upload, 'LEDGER', build / 'uploaded_hashes.json')
    return tmp_path, reg


class FakeGitHub:
    """An in-memory release store standing in for release_upload._request;
    every request is recorded."""

    def __init__(self, tags):
        self.calls = []
        self.releases = {}
        self.next_id = 100
        for t in tags:
            self._new(t, draft=False)

    def _new(self, tag, draft):
        self.next_id += 1
        self.releases[self.next_id] = {'id': self.next_id, 'tag_name': tag,
                                       'draft': draft, 'assets': {}}
        return self.releases[self.next_id]

    def id_of(self, tag):
        return next(i for i, r in self.releases.items() if r['tag_name'] == tag)

    def __call__(self, url, method='GET', data=None, content_type=None,
                 retries=4):
        self.calls.append((method, url, data))
        base = f'{API}/repos/{REPO}/releases'
        if url.startswith(base + '?') and method == 'GET':
            if 'page=1' not in url.split('&'):
                return []
            return [{k: v for k, v in r.items() if k != 'assets'}
                    for r in self.releases.values()]
        if url == base and method == 'POST':
            body = json.loads(data)
            r = self._new(body['tag_name'], body['draft'])
            return {k: v for k, v in r.items() if k != 'assets'}
        if '/assets?per_page' in url and method == 'GET':
            rid = int(url.split('/releases/')[1].split('/')[0])
            if 'page=1' not in url.split('&'):
                return []
            return list(self.releases[rid]['assets'].values())
        if url.startswith(f'{UPLOADS}/repos/{REPO}/releases/') and method == 'POST':
            rid = int(url.split('/releases/')[1].split('/')[0])
            name = url.split('?name=')[1]
            self.next_id += 1
            self.releases[rid]['assets'][name] = {'id': self.next_id,
                                                  'name': name, 'size': len(data)}
            return {}
        if url.startswith(base + '/assets/') and method == 'DELETE':
            aid = int(url.rsplit('/', 1)[1])
            for r in self.releases.values():
                r['assets'] = {n: a for n, a in r['assets'].items()
                               if a['id'] != aid}
            return None
        if method == 'PATCH':
            rid = int(url.rsplit('/', 1)[1])
            self.releases[rid].update(json.loads(data))
            return {}
        raise AssertionError('unexpected request {} {}'.format(method, url))

    def only_on(self, tag):
        """Every recorded request is the release-list GET, the creation of
        `tag`'s release, or a request on that release's id, upload URL or
        one of its assets."""
        rid = self.id_of(tag)
        aids = {a['id'] for a in self.releases[rid]['assets'].values()}
        base = f'{API}/repos/{REPO}/releases'
        for method, url, data in self.calls:
            if url.startswith(base + '?per_page=') and method == 'GET':
                continue
            if url == base and method == 'POST':
                assert json.loads(data)['tag_name'] == tag
                continue
            if url.startswith((f'{base}/{rid}/', f'{base}/{rid}?')) or \
                    url == f'{base}/{rid}':
                continue
            if url.startswith(f'{UPLOADS}/repos/{REPO}/releases/{rid}/assets?'):
                continue
            if url.startswith(f'{base}/assets/') and method == 'DELETE':
                assert int(url.rsplit('/', 1)[1]) not in {
                    a['id'] for r in self.releases.values()
                    if r['id'] != rid for a in r['assets'].values()}
                continue
            raise AssertionError('request off {}: {} {}'.format(tag, method, url))


def _gh(monkeypatch, tags=('v2', 'v2-endf')):
    gh = FakeGitHub(tags)
    monkeypatch.setattr(release_upload, '_request', gh)
    return gh


def test_plan_of_the_revision_tag(layout, capsys):
    release_upload.plan('v2.1')
    out = capsys.readouterr().out.split('\n')
    assert out[0] == 'v2.1: 3 assets'
    assert [s.strip() for s in out[1:] if s] == \
           ['decay.db', 'SHA256SUMS', 'data_manifest.toml']
    _, plan = release_upload._plan(tag='v2.1', need_meta=False)
    tmp, _ = layout
    assert plan == [('v2.1', [tmp / 'build' / 'decay.db',
                              tmp / 'build' / 'v2.1' / 'SHA256SUMS',
                              tmp / 'data_manifest.toml'])]
    assert not (tmp / 'build' / 'v2.1').exists()     # plan writes nothing


def test_generation_tag_plan_holds_no_revised_file(layout):
    tmp, reg = layout
    _, plan = release_upload._plan(need_meta=False)
    tags = dict(plan)
    names = [p.name for p in tags['v2']]
    assert 'decay.db' not in names
    assert names == sorted(f for f in reg['files'] if f != 'decay.db') + \
        ['SHA256SUMS', 'data_manifest.toml']
    assert [p.name for p in tags['v2.1']] == \
           ['decay.db', 'SHA256SUMS', 'data_manifest.toml']
    held = [t for t, paths in plan if 'decay.db' in [p.name for p in paths]]
    assert held == ['v2.1']


def test_without_revisions_the_whole_generation_is_planned_and_uploaded(
        layout, monkeypatch):
    """With no revisions the plan is the one-generation layout (every whole
    file + build/SHA256SUMS + manifest on the generation tag, then each
    shard group on its tag), and upload/verify with no tag run on all of
    it, writing build/SHA256SUMS."""
    tmp, reg = layout
    reg = copy.deepcopy(reg)
    del reg['revisions']
    release_upload.REGISTRY.write_text(json.dumps(reg))
    gen, plan = release_upload._plan()
    expected = [('v2', [tmp / 'build' / f for f in sorted(reg['files'])]
                 + [tmp / 'build' / 'SHA256SUMS', tmp / 'data_manifest.toml'])]
    for lib, group in sorted(reg['shards'].items()):
        expected.append((_tag(group['base_url']),
                         [tmp / 'shards-v2' / lib / f
                          for f in sorted(group['files'])]))
    assert gen == 'v2'
    assert plan == expected
    gh = _gh(monkeypatch)
    release_upload.upload()
    assert (tmp / 'build' / 'SHA256SUMS').read_text() == ''.join(
        '{}  {}\n'.format(reg['files'][f], f) for f in sorted(reg['files']))
    assert release_upload.verify()
    assert {r['tag_name'] for r in gh.releases.values()} == {t for t, _ in plan}


@pytest.mark.parametrize('command', ['upload', 'verify', 'publish'])
@pytest.mark.parametrize('tag', [None, 'v2', 'v2-endf', 'v9'])
def test_commands_off_the_revision_tag_exit_before_any_request(
        layout, monkeypatch, command, tag):
    tmp, _ = layout
    sums = (tmp / 'build' / 'SHA256SUMS').read_bytes()
    ledger = (tmp / 'build' / 'uploaded_hashes.json').read_bytes()
    gh = _gh(monkeypatch)
    with pytest.raises(SystemExit) as e:
        getattr(release_upload, command)(tag)
    assert 'revision tag' in str(e.value) or 'no tag' in str(e.value)
    assert gh.calls == []
    assert (tmp / 'build' / 'SHA256SUMS').read_bytes() == sums
    assert (tmp / 'build' / 'uploaded_hashes.json').read_bytes() == ledger
    assert not (tmp / 'build' / 'v2.1').exists()


def test_upload_verify_publish_touch_only_the_revision_release(
        layout, monkeypatch):
    tmp, reg = layout
    sums = (tmp / 'build' / 'SHA256SUMS').read_bytes()
    ledger0 = json.loads((tmp / 'build' / 'uploaded_hashes.json').read_text())
    gh = _gh(monkeypatch)
    release_upload.upload('v2.1')
    rid = gh.id_of('v2.1')
    assert gh.releases[rid]['draft']
    assert set(gh.releases[rid]['assets']) == \
        {'decay.db', 'SHA256SUMS', 'data_manifest.toml'}
    assert release_upload.verify('v2.1')
    release_upload.publish('v2.1')
    assert not gh.releases[rid]['draft']
    gh.only_on('v2.1')
    for r in gh.releases.values():                   # v2 and v2-endf untouched
        if r['id'] != rid:
            assert r['assets'] == {} and not r['draft']

    # SHA256SUMS of the revision tag: its one file, under build/v2.1/
    rsums = (tmp / 'build' / 'v2.1' / 'SHA256SUMS').read_text()
    assert rsums == '{}  decay.db\n'.format(reg['files']['decay.db'])
    assert gh.releases[rid]['assets']['SHA256SUMS']['size'] == len(rsums)
    # the generation's SHA256SUMS is byte-identical
    assert (tmp / 'build' / 'SHA256SUMS').read_bytes() == sums
    # the ledger gains the revision tag and keeps every other entry
    ledger = json.loads((tmp / 'build' / 'uploaded_hashes.json').read_text())
    assert {k: v for k, v in ledger.items() if k != 'v2.1'} == ledger0
    assert ledger['v2.1']['decay.db'] == reg['files']['decay.db']


def test_upload_onto_an_existing_revision_draft(layout, monkeypatch):
    gh = _gh(monkeypatch, tags=('v2', 'v2-endf', 'v2.1'))
    rid = gh.id_of('v2.1')
    gh.releases[rid]['draft'] = True
    gh.releases[rid]['assets']['decay.db'] = {'id': 1, 'name': 'decay.db',
                                              'size': 3}
    release_upload.upload('v2.1')
    release_upload.publish('v2.1')
    gh.only_on('v2.1')
    assert not any(c[0] == 'POST' and c[1].endswith('/releases')
                   for c in gh.calls)
    assert gh.releases[rid]['assets']['decay.db']['size'] == len(
        b'whole file decay.db\n')


def test_verify_hashes_the_revised_file_only(layout, monkeypatch):
    tmp, _ = layout
    gh = _gh(monkeypatch, tags=('v2', 'v2-endf', 'v2.1'))
    (tmp / 'build' / 'v2.1').mkdir()
    (tmp / 'build' / 'v2.1' / 'SHA256SUMS').write_text('x\n')
    hashed = []
    real = release_upload._sha256

    def spy(p):
        hashed.append(p.name)
        return real(p)
    monkeypatch.setattr(release_upload, '_sha256', spy)
    assert not release_upload.verify('v2.1')         # remote assets absent
    assert hashed == ['decay.db']
    gh.only_on('v2.1')


def test_verify_catches_a_revised_file_off_the_registry(layout, monkeypatch,
                                                        capsys):
    tmp, _ = layout
    _gh(monkeypatch)
    release_upload.upload('v2.1')
    (tmp / 'build' / 'decay.db').write_bytes(b'another decay.db\n')
    assert not release_upload.verify('v2.1')
    assert '1 checked, 1 mismatched' in capsys.readouterr().out


def _decay_db(path, generation):
    con = sqlite3.connect(path)
    con.execute('CREATE TABLE _version (library TEXT, generation TEXT)')
    con.execute('INSERT INTO _version VALUES (?,?)', ('decay', generation))
    con.commit()
    return con


def test_check_decay_stamp_follows_the_curie_registry(tmp_path):
    assert validate.check_decay_stamp(_decay_db(tmp_path / 'a.db', 'v2.1'),
                                      CURIE_REGISTRY)
    assert not validate.check_decay_stamp(_decay_db(tmp_path / 'b.db', 'v2'),
                                          CURIE_REGISTRY)
    # the default registry is the sibling checkout's
    assert validate.check_decay_stamp(_decay_db(tmp_path / 'c.db', 'v2.1'))
    # without a revision the top-level tag is expected
    reg = json.loads(CURIE_REGISTRY.read_text())
    del reg['revisions']
    regp = tmp_path / 'reg.json'
    regp.write_text(json.dumps(reg))
    assert validate.check_decay_stamp(_decay_db(tmp_path / 'd.db', 'v2'), regp)
    assert not validate.check_decay_stamp(_decay_db(tmp_path / 'e.db', 'v2.1'),
                                          regp)
