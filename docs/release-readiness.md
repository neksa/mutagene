# Prepare a public release

The source changes require a new tested release before users can install them from PyPI.
The existing public package does not include this draft's new workflows.

## Validate the artifact

Build the source distribution and wheel with `uv build`.
Install the wheel in a clean environment outside the checkout.
Run `tests/packaging/smoke_wheel.py` with the web extra installed.
Run the full test suite and documentation build.
Review the benchmark provenance and manuscript claims independently.

## Publish deliberately

Select the release version and update package metadata consistently.
Confirm repository ownership, release permissions, and the PyPI project configuration.
Configure a PyPI trusted publisher before enabling any package-upload workflow.
The existing release workflow does not establish PyPI publication.
Create an archival software record and link the manuscript to that exact release.

## Separate research artifacts

The benchmark uses third-party reference signatures from a pinned upstream repository.
Review applicable data terms before redistributing the evidence bundle publicly.
Keep AlphaMissense tables in the external cache and preserve their attribution manifest.
The workflow examples remain local publication candidates until their own CI and release records exist.
