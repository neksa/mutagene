import logging
from pathlib import Path


class DiscoverMenu:
    def __init__(self, parser):
        parser.add_argument(
            "-i",
            "--infile",
            type=Path,
            help="Integer count matrix with channels in rows and tumors in columns",
        )
        parser.add_argument(
            "--reference-matrix",
            type=Path,
            help="Reference signatures with matching channel labels",
        )
        parser.add_argument(
            "-o",
            "--output-dir",
            type=Path,
            help="Directory for candidate profiles, reference matches, and held-out errors",
        )
        parser.add_argument("--withhold", help="Omit one reference in a separate recovery control")

    def callback(self, args):
        from mutagene.signatures.discovery import compare_candidate

        if not all((args.infile, args.reference_matrix, args.output_dir)):
            raise SystemExit("Supply --infile, --reference-matrix, and --output-dir")
        try:
            summary = compare_candidate(
                args.infile, args.reference_matrix, args.output_dir, args.withhold
            )
        except (OSError, RuntimeError, ValueError) as error:
            logging.getLogger(__name__).error(str(error))
            raise SystemExit(1) from error
        print(summary.to_string(index=False))
