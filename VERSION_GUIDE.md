# Version Management

MutaGene uses semantic versioning (MAJOR.MINOR.PATCH) with GitHub Actions for automated version bumping.

## How to Bump Version

### Automatic (via PR labels):

1. Create a PR with your changes
2. Add one of these labels to the PR:
   - `patch` - Bug fixes (2.0.0 → 2.0.1)
   - `minor` - New features (2.0.0 → 2.1.0)
   - `major` - Breaking changes (2.0.0 → 3.0.0)
3. Merge the PR
4. GitHub Actions will automatically:
   - Update version in `pyproject.toml` and `mutagene/version.py`
   - Create a commit and git tag
   - Create a GitHub Release

### Manual (if needed):

```bash
# 1. Update version in both files
vim pyproject.toml  # Change: version = "2.0.1"
vim mutagene/version.py  # Change: __version__ = '2.0.1'

# 2. Commit and tag
git add pyproject.toml mutagene/version.py
git commit -m "chore: bump version to 2.0.1"
git tag v2.0.1
git push --tags
```

## For v2 Major Release

The v2-refactor branch PR should use the `major` label to bump from 0.9.2.0 → 2.0.0.

After the v2 release, stay in 2.x territory:
- Use `patch` for bug fixes (2.0.0 → 2.0.1)
- Use `minor` for new features (2.0.0 → 2.1.0)
- Reserve `major` for next breaking changes (2.x → 3.0.0)
