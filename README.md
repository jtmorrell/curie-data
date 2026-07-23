# curie-data

Versioned nuclear-data releases and data tooling for [curie](https://github.com/jtmorrell/curie),
a Python toolkit for the analysis of experimental nuclear data.

## Downloading the data files

curie downloads these files automatically on first use — no manual step is
needed on a machine with ordinary network access. If you are here to fetch
them by hand (for example on a network whose proxy blocks automatic
downloads), the current data generation's databases are:

| File | Library | Size |
|---|---|---|
| [`decay.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/decay.db) | Decay data: half-lives, emission energies and intensities, branching ratios (ENSDF via IAEA LiveChart, NUBASE2020/AME2020) | 117 MB |
| [`endf.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/endf.db) | ENDF/B-VIII.1 neutron cross sections | 842 MB |
| [`tendl.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/tendl.db) | TENDL-2025 neutron cross sections, exclusive channels | 37 MB |
| [`tendl_n_rp.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/tendl_n_rp.db) | TENDL-2025 neutron residual-product cross sections | 55 MB |
| [`tendl_p_rp.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/tendl_p_rp.db) | TENDL-2025 proton residual-product cross sections | 40 MB |
| [`tendl_d_rp.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/tendl_d_rp.db) | TENDL-2025 deuteron residual-product cross sections | 37 MB |
| [`tendl_a_rp.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/tendl_a_rp.db) | TENDL-2025 alpha residual-product cross sections | 32 MB |
| [`IRDFF.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/IRDFF.db) | IRDFF-II neutron dosimetry cross sections | 38 MB |
| [`iaea_monitors.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/iaea_monitors.db) | IAEA charged-particle monitor reactions (2025) | 10 MB |
| [`ziegler.db`](https://github.com/jtmorrell/curie-data/releases/download/v2/ziegler.db) | Stopping powers and atomic data | 0.2 MB |

Place the files you need in curie's data directory — `%LOCALAPPDATA%\curie` on
Windows (almost always `C:\Users\<username>\AppData\Local\curie`),
`~/.local/share/curie` on Linux, `~/Library/Application Support/curie` on
macOS — and curie uses them as-is. Downloads can be verified against
[`SHA256SUMS`](https://github.com/jtmorrell/curie-data/releases/download/v2/SHA256SUMS)
from the same release. The full manual-installation recipe (including
installing curie itself without PyPI access) is in the
[curie installation guide](https://jtmorrell.github.io/curie/quickinstall.html).

## Release layout

This repository's **Releases are releases of data, not software**: each tag is an
immutable-by-convention snapshot of the nuclear databases that a given range of curie
versions consumes. The curie package itself ships no data — it carries a registry
(`curie/data_registry.json` in the curie repo) mapping every file here to its SHA256 and
download URL, fetches files on first use, and verifies every download against the
registry checksum.

One data generation spans several release tags, because GitHub caps a release at 1000
assets and the sharded libraries carry hundreds of files each. The current generation
(`v2`) comprises:

| Tag | Contents |
|---|---|
| `v2` | The 10 whole database files + `SHA256SUMS` + `data_manifest.toml` |
| `v2-endf` | ENDF/B-VIII.1 per-target shards (`endf_<TABLE>.db`, one table each) |
| `v2-tendl` | TENDL-2025 exclusive-channel shards |
| `v2-tendl-n-rp` | TENDL-2025 neutron residual-product shards |
| `v2-tendl-p-rp` | TENDL-2025 proton residual-product shards |
| `v2-tendl-d-rp` | TENDL-2025 deuteron residual-product shards |
| `v2-tendl-a-rp` | TENDL-2025 alpha residual-product shards |

The `v1`/`v1-<library>` tags are the previous generation (ENDF/B-VII.1, TENDL-2015),
kept for the curie versions pinned to them.

The whole files and the shards are two representations of the same libraries: the shard
tooling pins both to the same source file in the registry, so they cannot silently
diverge. curie's lazy path fetches individual shards (a reaction lookup costs a fraction
of a MB); `ci.download()` fetches whole files for offline preparation.

Databases that are not sharded (decay, ziegler, IRDFF, iaea_monitors) are fetched whole:
decay.db is traversed globally by decay-chain calculations, IRDFF stores per-reaction
rather than per-target tables, and the others are small.

## Versioning and compatibility

- Every released curie version pins **exactly one** data generation through its
  registry: curie 0.1.0–0.2.x → `v1`; curie 0.3.0+ → `v2`. Old tags are never modified
  or deleted, so old curie versions keep working indefinitely.
- A new data generation is published as a new set of tags (`v3`, `v3-endf`, ...)
  alongside the old ones, and ships to users via a new curie release carrying the
  updated registry.
- A library whose shard count approaches the per-release asset cap splits by Z-range
  into multiple tags (`v3-tendl-1`, `v3-tendl-2`, ...); the registry's per-group
  base_url makes this transparent to curie.
- Historical note: curie 0.0.35–0.0.37 fetch whole files from the frozen `data-v1`
  release on the curie repository itself; that release stays in place for them.

## Tooling

The full pipeline that produces a data generation lives here:

- `scripts/fetch/` downloads the upstream evaluations (ENDF/B, TENDL, IAEA
  LiveChart, IRDFF-II, IAEA monitors), pinned in `sources_manifest.toml`.
- `scripts/build/` builds each database from the fetched sources. The TENDL
  libraries carry a curated target set — naturally occurring or half-life of at
  least one year, ground states and isomeric targets alike — defined by
  `scripts/build/tendl_targets.txt`, which is the single filter shared by fetch,
  build, and validation.
- `scripts/shard_library.py` splits a per-target library into one small sqlite
  file per table and records per-shard SHA256s in the curie registry (default:
  the sibling `../curie` checkout). Run it against the exact database file being
  published — it also pins the whole-file hash so the two representations stay
  coupled.
- `scripts/validate.py` runs the build-time checks — the schema contract
  curie's readers require, plus a coverage diff against the previous
  generation — any failure of which blocks release; `scripts/release_upload.py`
  publishes the release assets and verifies them against the registry
  checksums.

`data_manifest.toml` records the provenance (source library, size, SHA256) of the
current data generation.
