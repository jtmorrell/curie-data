"""Split a cross-section library into per-target shard files for a curie data release.

The ENDF and TENDL libraries store one table per target (e.g. CD_115, all
outgoing channels or residual products as columns) plus the small
all_reactions search index. Downloading a whole library to look up one target
is wasteful, so the data release carries one small sqlite file per table;
curie fetches only the shards a session actually touches and reassembles them
into its local data-directory database.

Usage (from a checkout of curie-data, with curie checked out alongside):
    python scripts/shard_library.py <library> [--source ../curie/curie/data/<library>.db]
                                              [--out shards/<library>]
                                              [--registry ../curie/curie/data_registry.json]

where <library> is one of: endf, tendl, tendl_n_rp, tendl_p_rp, tendl_d_rp, tendl_a_rp.

Writes <out>/<library>_<TABLE>.db, one per table, and updates the library's
shard group in the curie registry (per-shard SHA256s; the group's base_url —
its own release tag on this repository, so each release stays well under the
1000-assets-per-release limit — is preserved, or seeded for a new group). The
whole-file hash in the registry is pinned to the same source file, so the
eager and shard-assembled representations can never silently diverge. Table
names are validated as [A-Za-z0-9_] so the table<->filename mapping is
trivially invertible.
"""

import argparse
import hashlib
import json
import os
import re
import sqlite3


def sha256(path):
	h = hashlib.sha256()
	with open(path, 'rb') as f:
		for chunk in iter(lambda: f.read(1 << 20), b''):
			h.update(chunk)
	return h.hexdigest()


def shard(library, source, out_dir):
	os.makedirs(out_dir, exist_ok=True)
	# stale shards from a previous run would be uploaded as assets the
	# registry never references; match the full shard-name shape so that
	# e.g. 'tendl' cleanup cannot eat 'tendl_n_rp' shards in a shared out dir
	stale = re.compile(r'^{}_([A-Z][A-Za-z0-9_]*|all_reactions|_version)\.db$'.format(re.escape(library)))
	for old in os.listdir(out_dir):
		if stale.match(old):
			os.remove(os.path.join(out_dir, old))
	src = sqlite3.connect(source)
	# SQLite-internal tables (sqlite_sequence etc.) and schemaless rows are
	# not shardable content. Underscore tables are metadata: _version ships
	# as a shard (the assembled local database needs the generation stamp
	# for the connection-time skew check), build-internal ones (_build_meta)
	# live in the whole file only.
	schemas = {n: s for n, s in src.execute("SELECT name, sql FROM sqlite_master WHERE type='table'")
			   if s and not n.startswith('sqlite_')
			   and (not n.startswith('_') or n == '_version')}
	src.close()

	entries = {}
	for n, tb in enumerate(sorted(schemas)):
		if not re.fullmatch(r'[A-Za-z0-9_]+', tb):
			raise ValueError('table name {!r} is not filename-safe'.format(tb))
		fnm = '{}_{}.db'.format(library, tb)
		path = os.path.join(out_dir, fnm)
		con = sqlite3.connect(path)
		con.execute("ATTACH DATABASE ? AS src", (source,))
		con.execute(schemas[tb])
		con.execute('INSERT INTO "{0}" SELECT * FROM src."{0}"'.format(tb))
		con.commit()
		con.execute("DETACH DATABASE src")
		con.execute("VACUUM")
		con.close()
		entries[fnm] = sha256(path)
		if (n + 1) % 50 == 0 or n + 1 == len(schemas):
			print('{}/{} shards written'.format(n + 1, len(schemas)))
	return entries


def main():
	ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	ap.add_argument('library', choices=['endf', 'tendl', 'tendl_n_rp', 'tendl_p_rp', 'tendl_d_rp', 'tendl_a_rp'])
	ap.add_argument('--source', default=None)
	ap.add_argument('--out', default=None)
	ap.add_argument('--registry', default=os.path.join(root, '..', 'curie', 'curie', 'data_registry.json'))
	ap.add_argument('--release-prefix', default=None,
		help="release-tag prefix for this data generation (e.g. 'v2'): points the "
		     "group's base_url at <prefix>-<library>; omitted, an existing group's "
		     "base_url is preserved and a new one seeds at v1")
	args = ap.parse_args()
	source = args.source or os.path.join(root, '..', 'curie', 'curie', 'data', args.library + '.db')
	out = args.out or os.path.join(root, 'shards', args.library)

	registry_path = args.registry
	if not os.path.exists(registry_path):
		raise SystemExit('{} not found: pass --registry pointing at a curie checkout (the registry carries the base_url and whole-file hashes).'.format(registry_path))
	with open(registry_path) as f:
		registry = json.load(f)

	# the whole-file entry and the shards are two representations of the same
	# library; pin both to this source so they can never silently diverge
	source_sha = sha256(source)
	recorded = registry['files'].get(args.library + '.db')
	if recorded and recorded != source_sha:
		print('NOTE: whole-file {}.db hash updated ({}... -> {}...); publish the new {}.db alongside the shards.'.format(
			args.library, recorded[:12], source_sha[:12], args.library))
	registry['files'][args.library + '.db'] = source_sha

	entries = shard(args.library, source, out)
	group = registry.setdefault('shards', {}).setdefault(args.library, {})
	if args.release_prefix:
		url = 'https://github.com/jtmorrell/curie-data/releases/download/{}-{}/'.format(
			args.release_prefix, args.library.replace('_', '-'))
		if group.get('base_url') != url:
			group['base_url'] = url
			print('NOTE: shard group base_url set to {} - create that release tag before shipping.'.format(url))
	elif 'base_url' not in group:
		group['base_url'] = 'https://github.com/jtmorrell/curie-data/releases/download/v1-{}/'.format(args.library.replace('_', '-'))
		print('NOTE: new shard group; base_url set to {} - create that release tag before shipping.'.format(group['base_url']))
	group['files'] = entries
	with open(registry_path, 'w') as f:
		json.dump(registry, f, indent=1, sort_keys=True)
	print('registry updated: {} ({} {} shards)'.format(registry_path, len(entries), args.library))


if __name__ == '__main__':
	main()
