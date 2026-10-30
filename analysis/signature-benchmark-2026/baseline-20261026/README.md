# Preserve the pre-repair evidence

This directory preserves source and documentation snapshots from before the likelihood repairs.
The original manifest copies record the paths and hashes at that time.
`relocations.json` maps historical manifest entries to these snapshots.
Relocation changes paths, not artifact contents or recorded hashes.

Historical benchmark outputs remain in their original result directories.
The verification scripts resolve historical source hashes through `baseline_sources.py`.
The new fresh-panel comparison imports the archived `identify.py` for pre-repair MLE.
That comparison explicitly requests MLE and does not invoke archived sparse fitting.

To reproduce the original sparse search or likelihood audit, use an isolated checkout with the archived fitter sources.
Do not restore those sources over the repaired production files during a benchmark run.
The repaired benchmark records its own source hashes before fitting starts.
