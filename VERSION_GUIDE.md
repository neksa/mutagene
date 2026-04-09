# Version Management

MutaGene uses semantic versioning (MAJOR.MINOR.PATCH).

Version is defined in `pyproject.toml` and `mutagene/version.py`.

## Bumping via PR labels (preferred)

Add a label to your PR before merging:
- `patch` — bug fixes (1.0.0 -> 1.0.1)
- `minor` — new features (1.0.0 -> 1.1.0)
- `major` — breaking changes (1.0.0 -> 2.0.0)

GitHub Actions updates both files, creates a tag and release on merge.
