# Prepare a release

Choose the version according to the compatibility rules in [CONTRIBUTING.md](CONTRIBUTING.md).
Change the version in `pyproject.toml` and `mutagene/version.py`.
Run `uv lock` and add release notes to the matching section in `CHANGELOG.md`.
Add the corresponding `major`, `minor`, or `patch` label to the release pull request.

After merge, GitHub Actions tags the merge commit and creates the GitHub release.
PyPI publication requires a separate upload.
