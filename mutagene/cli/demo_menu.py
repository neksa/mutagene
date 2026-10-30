"""Run an offline demonstration with bundled synthetic mutations."""

import json
import logging
import shutil
from pathlib import Path

from mutagene.io.context_window import read_mutations
from mutagene.io.decomposition import write_decomposition
from mutagene.io.matrix import write_count_matrix
from mutagene.io.profile import read_signatures
from mutagene.profiles.profile import get_multisample_mutational_profile
from mutagene.signatures.identify import decompose_mutational_profile_counts
from mutagene.version import __version__

logger = logging.getLogger(__name__)


class DemoMenu:
    def __init__(self, parser):
        parser.description = "Run a small synthetic demonstration without network access."
        parser.add_argument("--output-dir", type=Path, default=Path("mutagene-demo"))

    def callback(self, args):
        destination = args.output_dir
        if destination.exists() and any(destination.iterdir()):
            raise SystemExit("Choose an empty output directory for the demo")
        destination.mkdir(parents=True, exist_ok=True)
        source = Path(__file__).resolve().parents[1] / "data" / "demo"
        for name in ("sample.maf", "reference.2bit"):
            shutil.copyfile(source / name, destination / name)
        with (destination / "sample.maf").open() as stream:
            mutations, _, stats = read_mutations(
                "MAF", stream, str(destination / "reference.2bit"), window_size=1
            )
        profiles = get_multisample_mutational_profile(mutations, counts=True)
        with (destination / "matrix.tsv").open("w") as stream:
            write_count_matrix(stream, profiles)
        signatures = read_signatures("COSMICv2")
        results = {}
        for sample, profile in profiles.items():
            _, _, results[sample] = decompose_mutational_profile_counts(
                profile, signatures, "MLE", others_threshold=0, enable_dummy=True
            )
        with (destination / "signatures.tsv").open("w") as stream:
            write_decomposition(stream, results, signatures[1])
        summary = {
            "mutagene_version": __version__,
            "synthetic": True,
            "samples": len(profiles),
            "mutations": stats["loaded"],
            "channels": 96,
            "interpretation": "This tiny fixture demonstrates execution, not biological inference.",
        }
        (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(f"Demo complete: {destination.resolve()}")
        print("Expected result: 3 samples, 24 mutations, 96 channels.")
