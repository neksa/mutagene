# Preserve the likelihood reassessment

These snapshots retain files that change during the background-selection assessment.
The original manifest records their original paths and hashes.
The relocated manifest retains the same hashes at archived paths.
`relocations.json` maps those original paths to the snapshots.
Unchanged artifacts remain at their original locations.

Restore these snapshots in an isolated checkout before repeating the preceding benchmark.
The preceding reassessment uses full-model MLE as its interface default.
The background report records the later method-selection decision.
