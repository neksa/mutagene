# Changelog

Notable changes to MutaGene. Dates are release dates; the format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and versions follow
[Semantic Versioning](https://semver.org/).

## [1.1.1] - 2026-09-03

### Fixed

- The published documentation at
  [mutagene.readthedocs.io](https://mutagene.readthedocs.io) had not been
  rebuilt: it described 1.0.0, returned 404 for pages added since, and still
  carried the incorrect claim that `rank` uses a pan-cancer cohort by default.
  Read the Docs has required a configuration file since 2023 and there was none,
  so no build ran. The documented version is now read from the package rather
  than repeated, and the build runs in CI so a failure is visible in the pull
  request rather than only on Read the Docs.

[1.1.1]: https://github.com/neksa/mutagene/compare/v1.1.0...v1.1.1

## [1.1.0] - 2026-08-31

The first release since the web interface landed. Most of this is correctness
work on the mutation parsers and the ranking statistics, where the failures were
silent: analyses completed and returned plausible numbers that were wrong.

### Behaviour changes

These change what existing commands produce. None is a bug; each is a decision.

- **Variants the caller rejected are no longer counted.** A `FILTER` column is
  now honoured, so only `PASS` variants (and files that record no filters) are
  read. Mutation counts drop for inputs that carry the column. Use
  `--keep-filtered` for the previous behaviour, and `--filter-column NAME` for
  converters that write the verdict elsewhere.
- **`rank` output gained a `transcript` column and a `#` provenance header.**
  Readers of the TSV must skip comment lines: `pd.read_csv(..., comment="#")`.
  A `<outfile>.provenance.json` sidecar is written alongside.
- **A run that reads no mutations exits non-zero.** `profile`, `rank`, `motif`
  and `signature` previously wrote an empty output file and reported success.

### Added

- `--params-out` / `--params-in` on every subcommand, recording a run's
  parameters as JSON and replaying them, with command-line arguments taking
  precedence.
- `--keep-filtered` and `--filter-column` on the analysis subcommands.
- Genome-assembly mismatch detection, reported by `profile -v` and shown as a
  warning in the web interface.
- `serve` documentation, and a shared page for the arguments common to the
  analysis subcommands.

### Fixed

- MAF files carrying `Tumor_Seq_Allele2` without `Tumor_Seq_Allele1` were
  rejected by `rank`, `motif` and the web interface, though `profile` read them.
- One sample with no usable genomic context discarded every remaining sample in
  a multi-sample file, reporting the whole file as empty.
- The 5' and 3' bases were complemented in the wrong order for reference alleles
  reported on the reverse strand, putting mutations in the wrong channel of the
  96-channel profile.
- The genome-assembly mismatch warning could never fire, so analyses run against
  the wrong assembly completed silently.
- `rank` labelled a mutation present in every cohort sample "Undefined" instead
  of "Driver", because the sample's pseudocount pushed the observed count past
  the cohort size.
- Signature decomposition charged nothing for a mixture that assigns zero
  probability to a channel that was actually observed.
- `rank` chose the first transcript it encountered for each gene, so results
  depended on the order of rows in the MAF. The longest transcript is used now,
  with a deterministic fallback.
- A blank annotation column overwrote a populated one, silently dropping rows
  that had a usable annotation.
- The web interface offered VCF uploads and then read them as MAF.
- Re-running an analysis destroyed the previous results before the new ones were
  stored, so a failed re-run lost both.
- Unreadable input raised tracebacks out of the parsers rather than reporting
  what was wrong; truncated uploads named neither the file nor a remedy.
- Profile files accepted negative, NaN and infinite values, and silently kept
  the last of a duplicated channel.
- Resumed downloads requested one byte past the end of the file and never closed
  the response.
- `mutagene profile` had no error handling at all.
- The TCGI reader required `CHROM` while documenting `CHR`.

### Removed

- `mutagene/io/ensembl.py`, unreachable and carrying its own copy of the
  reverse-strand bug fixed above.

### Known issues

- `mutagene benchmark` cannot run in any mode ([#114]).
- `mutagene/reports/nci60.py` cannot execute ([#115]).

[1.1.0]: https://github.com/neksa/mutagene/compare/v1.0.0...v1.1.0
[#114]: https://github.com/neksa/mutagene/issues/114
[#115]: https://github.com/neksa/mutagene/issues/115

## [1.0.0] - 2025-11-27

Web interface, project modernisation and dependency upgrades. See the
[release notes](https://github.com/neksa/mutagene/releases/tag/v1.0.0).

[1.0.0]: https://github.com/neksa/mutagene/releases/tag/v1.0.0
