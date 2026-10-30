"""Verify the files included in the source-and-summary review snapshot."""

import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "analysis/publication-manifest.json").read_text())
    failures = []
    for name, record in manifest["files"].items():
        path = root / name
        if not path.is_file():
            failures.append(f"Missing file: {name}")
            continue
        content = path.read_bytes()
        if (
            len(content) != record["bytes"]
            or hashlib.sha256(content).hexdigest() != record["sha256"]
        ):
            failures.append(f"Content differs: {name}")
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"Verified {len(manifest['files'])} source-and-summary files")
    print("Full experiment verification requires the separately listed local artifacts")


if __name__ == "__main__":
    main()
