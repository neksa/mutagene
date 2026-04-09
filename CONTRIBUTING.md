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
