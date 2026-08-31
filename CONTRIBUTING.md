# Contributing

## Setup

1. Fork and clone: https://github.com/neksa/mutagene
2. Install [uv](https://docs.astral.sh/uv/getting-started/installation/)
3. `uv sync --dev` (add `--extra web` for webapp development)

## Workflow

1. Create a branch: `git checkout -b feature/my-change`
2. Make changes, add tests for new functionality
3. Verify:
   ```bash
   uv run pytest
   uv run ruff check mutagene tests
   uv run black --check mutagene tests
   ```
4. Submit a pull request (release notes go in the PR description)

## Releases

Releases are cut by merging a labelled pull request. Label every PR with one of:

| Label   | When |
|---------|------|
| `patch` | Bug fixes, documentation, internal changes. Nothing a user has to react to. |
| `minor` | New arguments or commands, and behaviour changes a user should read about before upgrading. |
| `major` | Removing or renaming something people depend on, or changing a default such that existing invocations produce materially different results. Deliberate, and not a substitute for a well-explained `minor`. |

An unlabelled PR merges without cutting a release, which is the right outcome
for work that is not worth releasing on its own.

Before merging a labelled PR, add the matching section to `CHANGELOG.md`:

```markdown
## [1.2.0] - 2026-09-14

### Behaviour changes
### Added
### Fixed
### Removed
```

The version-bump workflow refuses to release without that section, and the
GitHub release notes are taken from it. Put anything that changes what an
existing command produces under **Behaviour changes**, whether or not it is a
bug fix — that is the part people need before they upgrade.
