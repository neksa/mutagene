import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
from tqdm import tqdm

from mutagene.cli.genome_resolver import GenomeNotFoundError, resolve_genome
from mutagene.io.decomposition import write_decomposition
from mutagene.io.matrix import read_count_matrix, read_reference_matrix
from mutagene.io.profile import read_signatures
from mutagene.io.text import detect_mutation_format, open_text_input
from mutagene.io.variant_filter import add_filter_argument
from mutagene.profiles.profile import (
    generate_resampled_profiles,
    read_sample_profiles,
)
from mutagene.signatures.identify import IDENTIFY_MIN_FUNCTIONS, decompose_mutational_profile_counts
from mutagene.signatures.likelihood import DEFAULT_BACKGROUND_SELECTION, DEFAULT_METHOD
from mutagene.version import __version__

logger = logging.getLogger(__name__)

genome_error_message = """requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.soe.ucsc.edu/downloads.html for more
                          Use mutagene fetch to download genome assemblies"""


class SignatureMenu:
    def __init__(self, parser):
        parser.description = """
Identify the activity of known mutational signatures in samples.

Mutational signatures describe the distribution of mutations across channels. This command decomposes the mutational profile of
each sample into a non-negative linear combination of the signatures in the
chosen collection and reports the exposure of each signature.

For SBS1536 or ID83, supply a matching labelled reference matrix. ID83 also requires classified counts.

SBS96 signature collections available with --signatures (-s):

  MGA, 5        MutaGene-5, 5 signatures derived by MutaGene
  MGB, 10       MutaGene-10, 10 signatures derived by MutaGene
  COSMICv2, 30  COSMIC v2, 30 signatures
  COSMICv3, 49  COSMIC v3, 67 single-base-substitution signatures (SBS1-SBS85)
  KUCAB, 53     54 profiles of environmental mutagen exposure, including a
                Control profile

The numeric aliases are historical names and do not always match the number of
signatures in the collection.

Read more about MutaGene signatures at
https://www.ncbi.nlm.nih.gov/research/mutagene/signatures and about COSMIC
signatures at https://cancer.sanger.ac.uk/cosmic/signatures/

The decomposition is an optimization problem; --method (-m) selects the
function that is minimized. Confidence intervals for the exposures can be
obtained by bootstrap resampling with --bootstrap.
"""
        parser.epilog = """
Examples:
# identify COSMIC v2 signatures in sample1.maf using hg19
mutagene signature -i sample1.maf -g hg19 -s COSMICv2

# identify MutaGene-5 signatures and write the result to a file
mutagene signature -i sample1.maf -g hg19 -s 5 -o exposures.tsv

# add bootstrap confidence intervals from 1000 replicates
mutagene signature -i sample1.maf -g hg19 -s COSMICv3 --bootstrap -br 1000
        """

        required_group = parser.add_argument_group("Required arguments")
        required_group.add_argument(
            "--infile", "-i", help="Input file in VCF or MAF format", type=open_text_input
        )
        required_group.add_argument(
            "--genome",
            "-g",
            help="Location of genome assembly file in 2bit format",
            type=str,
            default="hg19",
        )
        catalogs = required_group.add_mutually_exclusive_group()
        catalogs.add_argument(
            "--signatures",
            "-s",
            choices=["5", "10", "30", "49", "53", "MGA", "MGB", "COSMICv2", "COSMICv3", "KUCAB"],
            help="Collection of signatures to use",
            type=str,
            default=None,
        )

        catalogs.add_argument(
            "--reference-matrix",
            type=open_text_input,
            metavar="FILE",
            help="Labelled channel-by-signature reference TSV",
        )
        optional_group = parser.add_argument_group("Optional arguments")
        optional_group.add_argument(
            "--context",
            choices=["SBS96", "SBS1536", "ID83"],
            default="SBS96",
            help="Mutation channels for fitting",
        )
        optional_group.add_argument(
            "--input-format",
            "-f",
            help="Input format: AUTO, MAF, VCF, TCGI, or a labelled count MATRIX",
            type=str.upper,
            choices=["AUTO", "MAF", "VCF", "TCGI", "MATRIX"],
            default="AUTO",
        )
        optional_group.add_argument(
            "--outfile",
            "-o",
            nargs="?",
            type=argparse.FileType("w"),
            default=sys.stdout,
            help="Name of output file, will be generated in TSV format",
        )

        # for backwards compatibility with 0.8.X add a hidden action that would just take anything as a valid input
        optional_group.add_argument("action", nargs="?", metavar="")

        add_filter_argument(optional_group)

        advanced_group = parser.add_argument_group("Advanced arguments")
        advanced_group.add_argument(
            "--method",
            "-m",
            help="Fitting objective or subset selection (default: sparse-bic)",
            type=str,
            default=DEFAULT_METHOD,
        )
        advanced_group.add_argument(
            "--background-selection",
            choices=["fixed", "joint", "block"],
            default=DEFAULT_BACKGROUND_SELECTION,
            help="Sparse background selection: fixed, individual components, or complete block",
        )
        advanced_group.add_argument(
            "--no-unexplained-variance",
            "-U",
            help="Do not account for unexplained variance (non-context dependent mutational processes and unknown signatures)",
            action="store_false",
        )
        advanced_group.add_argument(
            "--mutations-threshold",
            "-t",
            help="Only report signatures with mutations above the threshold",
            type=int,
            default=0,
        )
        advanced_group.add_argument(
            "--keep-only",
            "-k",
            help="Keep only the signatures in the list, separated by commas e.g. 1,3,5",
            type=str,
            default=None,
        )

        bootstrap_group = parser.add_argument_group("Bootstrap-specific arguments")
        bootstrap_group.add_argument(
            "--bootstrap",
            "-b",
            help="Use the bootstrap to calculate confidence intervals",
            action="store_true",
        )
        bootstrap_group.add_argument(
            "--bootstrap-replicates",
            "-br",
            help="Number of bootstrap replicates",
            type=int,
            default=100,
        )
        bootstrap_group.add_argument(
            "--bootstrap-confidence-level", "-bcl", help="Confidence level", type=int, default=95
        )
        bootstrap_group.add_argument(
            "--bootstrap-method",
            "-bm",
            help="Bootstrap method (t: t-distribution, p: percentile)",
            type=str,
            choices=["t", "p"],
            default="p",
        )
        optional_group.add_argument(
            "--diagnostics-out",
            type=str,
            help="Write sample counts and fitted null weights as JSON",
        )
        bootstrap_group.add_argument(
            "--seed", type=int, help="Seed bootstrap sampling for reproducible runs"
        )
        self.parser = parser

    def identify(self, args):
        if not args.infile:
            logger.warning(
                "Provide input file in VCF or MAF format (-i) and a corresponding genome assembly (-g)"
            )
            sys.exit(1)
        if args.context != "SBS96" and not args.reference_matrix:
            self.parser.error(f"{args.context} fitting requires --reference-matrix")
        if args.context == "ID83" and args.input_format != "MATRIX":
            self.parser.error("ID83 fitting requires a classified indel count matrix (-f MATRIX)")
        if args.input_format != "MATRIX":
            try:
                args.genome = resolve_genome(args.genome)
            except GenomeNotFoundError as e:
                logger.error(str(e))
                sys.exit(1)
        if args.method.lower() not in IDENTIFY_MIN_FUNCTIONS:
            logger.warning("Unknown method provided")
            sys.exit(1)

        if args.bootstrap_replicates < 10:
            logger.error("Specify at least 10 bootstrap replicates")
            sys.exit(1)

        if not 70 <= args.bootstrap_confidence_level < 99:
            logger.error("Specify a confidence level from 70% to less than 99%")
            sys.exit(1)

        only = None
        if args.keep_only is not None:
            only = [name.strip() for name in args.keep_only.split(",") if name.strip()]
            if len(only) < 1:
                logger.error("List of signatures to keep for the analysis is empty")
                sys.exit(1)
            logger.warning(
                "We will only analyze signatures in this list: {}".format(", ".join(only))
            )

        try:
            if args.reference_matrix is not None:
                W, signature_names = read_reference_matrix(args.reference_matrix, args.context)
                signature_set = args.reference_matrix.name
                if only:
                    if not set(only).issubset(signature_names) or len(only) != len(set(only)):
                        raise ValueError(
                            "The signature selection contains unknown or duplicate identifiers"
                        )
                    W = W[:, [signature_names.index(name) for name in only]]
                    signature_names = only
            else:
                signature_set = args.signatures or "COSMICv2"
                W, signature_names = read_signatures(signature_set, only=only)
                if not signature_names or (only and len(only) != len(signature_names)):
                    raise ValueError("The signature selection contains unknown identifiers")
            if args.input_format == "MATRIX":
                samples_profiles = read_count_matrix(args.infile, args.context)
            else:
                if args.input_format == "AUTO":
                    args.input_format = detect_mutation_format(args.infile)
                samples_profiles = read_sample_profiles(
                    args.input_format,
                    args.infile,
                    args.genome,
                    args.context,
                    keep_filtered=args.keep_filtered,
                    filter_column=args.filter_column,
                )
                if not samples_profiles:
                    raise ValueError("No mutations could be read for signature fitting")
        except (ValueError, OSError, EOFError) as error:
            logger.error(f"Cannot read {args.input_format} input: {error}")
            sys.exit(1)

        samples_results = {}
        diagnostics = {
            "mutagene_version": __version__,
            "signature_set": signature_set,
            "context": args.context,
            "method": args.method,
            "seed": args.seed,
            "samples": {},
        }
        rng = np.random.default_rng(args.seed)

        enable_background = args.context == "SBS96" and args.no_unexplained_variance
        for sample, profile in samples_profiles.items():
            optimizer = {}
            try:
                weights, _, results = decompose_mutational_profile_counts(
                    profile,
                    (W, signature_names),
                    args.method,
                    others_threshold=0.0,
                    enable_dummy=enable_background,
                    background_selection=args.background_selection,
                    diagnostics=optimizer,
                )
            except (ValueError, RuntimeError) as error:
                logger.error(f"Signature fitting fails for {sample}: {error}")
                sys.exit(1)
            samples_results[sample] = results
            diagnostics["samples"][sample] = {
                **optimizer,
                "mutations": int(np.sum(profile)),
                "known_signature_weight": float(np.sum(weights[: len(signature_names)])),
                "null_weight": float(np.sum(weights[len(signature_names) :])),
                "null_components": weights[len(signature_names) :].tolist(),
            }

        if not args.bootstrap:
            write_decomposition(
                args.outfile,
                samples_results,
                signature_names,
                mutations_threshold=args.mutations_threshold,
            )
        else:
            bootstrap_samples_results = {}
            for sample, profile in samples_profiles.items():
                bootstrap_results = []
                for resampled_profile in tqdm(
                    generate_resampled_profiles(profile, args.bootstrap_replicates, rng=rng),
                    total=args.bootstrap_replicates,
                ):
                    try:
                        _, _, results = decompose_mutational_profile_counts(
                            resampled_profile,
                            (W, signature_names),
                            args.method,
                            others_threshold=0.0,
                            enable_dummy=enable_background,
                            background_selection=args.background_selection,
                        )
                    except (ValueError, RuntimeError) as error:
                        logger.error(f"Bootstrap fitting fails for {sample}: {error}")
                        sys.exit(1)
                    bootstrap_results.append(results)
                bootstrap_samples_results[sample] = bootstrap_results

            write_decomposition(
                args.outfile,
                samples_results,
                signature_names,
                mutations_threshold=args.mutations_threshold,
                bootstrap_method=args.bootstrap_method,
                profile=samples_profiles,
                bootstrap_results=bootstrap_samples_results,
                bootstrap_level=args.bootstrap_confidence_level,
            )

        if args.diagnostics_out:
            Path(args.diagnostics_out).write_text(json.dumps(diagnostics, indent=2) + "\n")

    def callback(self, args):
        self.identify(args)
