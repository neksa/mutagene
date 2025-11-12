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

## During v2 Refactoring

Stay in 2.0.x territory for all refactoring work:
- Use `patch` for most changes
- Reserve `minor` for complete features
- Save `major` for post-refactor release
