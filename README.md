# curie-data

Versioned nuclear-data releases and data tooling for [curie](https://github.com/jtmorrell/curie),
a Python toolkit for the analysis of experimental nuclear data.

This repository's **Releases are releases of data, not software**: each tag is an
immutable-by-convention snapshot of the nuclear databases that a given range of curie
versions consumes. The curie package itself ships no data — it carries a registry
(`curie/data_registry.json` in the curie repo) mapping every file here to its SHA256 and
download URL, fetches files on first use, and verifies every download against the
registry checksum.

## Release layout

One data generation spans several release tags, because GitHub caps a release at 1000
assets and the sharded libraries carry ~420 files each:

| Tag | Contents |
|---|---|
| `v1` | The 9 whole database files + `SHA256SUMS` + `data_manifest.toml` |
| `v1-endf` | ENDF/B-VII.1 per-target shards (`endf_<TABLE>.db`, one table each) |
| `v1-tendl` | TENDL-2015 exclusive-reaction shards |
| `v1-tendl-n-rp` | TENDL-2015 neutron residual-product shards |
| `v1-tendl-p-rp` | TENDL-2015 proton residual-product shards |
| `v1-tendl-d-rp` | TENDL-2015 deuteron residual-product shards |

The whole files and the shards are two representations of the same libraries: the shard
tooling pins both to the same source file in the registry, so they cannot silently
diverge. curie's lazy path fetches individual shards (a reaction lookup costs a fraction
of a MB); `ci.download()` fetches whole files for offline preparation.

Databases that are not sharded (decay, ziegler, IRDFF, iaea_monitors) are fetched whole:
decay.db is traversed globally by decay-chain calculations, IRDFF stores per-reaction
rather than per-target tables, and the others are small.

## Versioning and compatibility

- Every released curie version pins **exactly one** data generation through its
  registry: curie 0.1.0 → `v1`. Old tags are never modified or deleted, so old curie
  versions keep working indefinitely.
- A new data generation (e.g. the ENDF/B-VIII.1 / TENDL-2023 rebuild) is published as a
  new set of tags (`v2`, `v2-endf`, ...) alongside the old ones, and ships to users via
  a new curie release carrying the updated registry.
- A library whose shard count approaches the per-release asset cap splits by Z-range
  into multiple tags (`v2-tendl-1`, `v2-tendl-2`, ...); the registry's per-group
  base_url makes this transparent to curie.
- Historical note: curie 0.0.35–0.0.37 fetch whole files from the frozen `data-v1`
  release on the curie repository itself; that release stays in place for them.

## Tooling

`scripts/shard_library.py` splits a per-target library (endf, tendl, tendl_n_rp,
tendl_p_rp, tendl_d_rp) into one small sqlite file per table and records per-shard
SHA256s in the curie registry (default: the sibling `../curie` checkout). Run it against
the exact database file being published — it also pins the whole-file hash so the two
representations stay coupled.

`data_manifest.toml` records the provenance (source library, size, SHA256) of the
current data generation.

The build pipeline that regenerates the databases themselves from the upstream
evaluations (ENDF/TENDL/ENSDF sources) will live here as it is developed.
