"""Analysis pipeline integration for webapp."""

import gzip
import json
import logging
from pathlib import Path
from typing import Any

from mutagene.io.context_stats import assembly_mismatch_warning
from mutagene.profiles.profile import calc_profile
from mutagene.version import __version__

from .genome_manager import GenomeManager

logger = logging.getLogger(__name__)

# A few hundred MB of mutation text is far beyond any realistic MAF; the cap
# stops a small tarball from expanding into a disk-filling extract.
MAX_EXTRACTED_BYTES = 2 * 1024 * 1024 * 1024

# Per-sample decomposition is linear but clustering builds an N x N distance
# matrix, so a file with many unique sample IDs is a memory-exhaustion vector.
MAX_CLUSTERING_SAMPLES = 1000


def extract_input_file(file_path: Path, output_dir: Path) -> Path:
    """Extract tar.gz archives to a stable location; return path to the mutation file.

    For non-archive files the original path is returned unchanged.
    Extracted files go into *output_dir* so they are cleaned up with the analysis.
    """
    import tarfile

    if str(file_path).endswith(".tar.gz") or str(file_path).endswith(".tgz"):

        def truncated(e):
            # Truncation surfaces from wherever the archive ran out -- opening
            # it, listing it, or streaming a member -- so every one of those has
            # to say the same actionable thing rather than leak an EOFError.
            return TruncatedInputError(
                f"{Path(file_path).name} could not be read as a tar.gz archive ({e}). "
                "The upload may have been interrupted; try uploading it again"
            )

        try:
            tar = tarfile.open(file_path, "r:gz")
        except (tarfile.ReadError, EOFError, gzip.BadGzipFile, OSError) as e:
            raise truncated(e) from None

        with tar:
            try:
                members = tar.getmembers()
            except (tarfile.ReadError, EOFError, gzip.BadGzipFile, OSError) as e:
                raise truncated(e) from None

            maf_candidates = [
                m
                for m in members
                if m.isfile() and ("mutation" in m.name.lower() or m.name.endswith(".maf"))
            ]
            if not maf_candidates:
                raise ValueError("No mutation file found in tar.gz archive")

            maf_member = maf_candidates[0]
            if maf_member.size > MAX_EXTRACTED_BYTES:
                raise ValueError(
                    f"Archived mutation file is too large "
                    f"({maf_member.size} bytes, limit {MAX_EXTRACTED_BYTES})"
                )
            logger.info(f"Extracting {maf_member.name} from tarball")

            extract_dir = output_dir / "extracted"
            extract_dir.mkdir(parents=True, exist_ok=True)

            # Never trust the member name: it may contain path separators, "..",
            # or be an absolute path. Stream the member out under a basename we
            # control instead of letting tarfile choose the destination.
            safe_name = Path(maf_member.name).name
            if not safe_name or safe_name in (".", ".."):
                raise ValueError(
                    f"Refusing to extract tar member with unsafe name: {maf_member.name!r}"
                )

            extracted_path = extract_dir / safe_name
            resolved = extracted_path.resolve()
            if not resolved.is_relative_to(extract_dir.resolve()):
                raise ValueError(
                    f"Tar member {maf_member.name} would extract outside target directory"
                )

            source = tar.extractfile(maf_member)
            if source is None:
                raise ValueError(f"Tar member {maf_member.name} is not a regular file")

            # Enforce the cap while copying rather than trusting the header size,
            # which an attacker controls independently of the actual stream.
            written = 0
            try:
                with source, open(extracted_path, "wb") as dest:
                    while chunk := source.read(1024 * 1024):
                        written += len(chunk)
                        if written > MAX_EXTRACTED_BYTES:
                            raise ValueError(
                                f"Archived mutation file exceeds {MAX_EXTRACTED_BYTES} bytes"
                            )
                        dest.write(chunk)
            except (tarfile.ReadError, EOFError, gzip.BadGzipFile) as e:
                extracted_path.unlink(missing_ok=True)
                raise truncated(e) from None
            except Exception:
                extracted_path.unlink(missing_ok=True)
                raise
            return extracted_path

    return file_path


class TruncatedInputError(ValueError):
    """An uploaded archive ended before it should have."""


def open_input_file(file_path: Path, mode: str = "rt"):
    """Open a mutation file, handling gzip transparently."""
    if str(file_path).endswith(".gz"):
        return _GuardedGzipFile(file_path, mode)
    return open(file_path, mode, encoding="utf-8")


class _GuardedGzipFile:
    """A gzip file whose truncation is reported in terms the uploader can act on.

    A part-uploaded .gz raises EOFError from wherever it happened to run out,
    which surfaced as "Compressed file ended before the end-of-stream marker was
    reached" with no indication of which file or what to do about it.
    """

    def __init__(self, file_path, mode):
        self._path = file_path
        self._fh = gzip.open(file_path, mode, encoding="utf-8")

    def _guard(self, call, *args):
        try:
            return call(*args)
        except (EOFError, gzip.BadGzipFile, OSError) as e:
            raise TruncatedInputError(
                f"{Path(self._path).name} is not a complete gzip file ({e}). "
                "The upload was probably interrupted; try uploading it again"
            ) from None

    def __iter__(self):
        # Iterating lazily, so truncation shows up part way through the read.
        iterator = iter(self._fh)
        while True:
            try:
                yield self._guard(next, iterator)
            except StopIteration:
                return

    def read(self, *args):
        return self._guard(self._fh.read, *args)

    def readlines(self, *args):
        return self._guard(self._fh.readlines, *args)

    def seek(self, *args):
        return self._fh.seek(*args)

    def close(self):
        self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def detect_input_format(file_path: Path) -> str:
    """Whether to read this upload as a MAF or a VCF.

    The profile step sniffs the format itself, but the per-sample read did not
    and asked for MAF unconditionally, so a VCF upload -- which the page offers
    and the extension allow-list accepts -- was profiled correctly and then died
    on "Chromosome is not defined in MAF file".
    """
    name = str(file_path).lower()
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    if name.endswith(".vcf"):
        return "VCF"
    if name.endswith(".maf"):
        return "MAF"

    # .txt and .tsv are allowed too, so fall back to the contents: a VCF opens
    # with ## metadata lines and a #CHROM header, neither of which a MAF has.
    try:
        with open_input_file(file_path, "rt") as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("##fileformat=VCF") or stripped.startswith("#CHROM"):
                    return "VCF"
                if not stripped.startswith("#"):
                    break
    except (OSError, ValueError):
        pass
    return "MAF"


def profile_channel_order() -> list[str]:
    """Return the canonical 96 mutation channel labels in signature-matrix row order.

    The signature matrices returned by ``read_signatures`` and the profiles written
    by ``calc_profile`` both follow ``get_profile_attributes_dict()`` order
    (5' base, 3' base, then mutation). Sorting the labels alphabetically yields a
    *different* permutation, so the two must never be zipped together by sorted key.
    """
    from mutagene.io.profile import get_profile_attributes_dict

    return [
        f"{a['context'][0]}[{a['mutation'][0]}>{a['mutation'][1]}]{a['context'][1]}"
        for a in get_profile_attributes_dict()
    ]


def profile_dict_to_array(profile_data: dict[str, int]):
    """Convert a ``{'A[C>T]G': count}`` dict into a 96-element array in signature order."""
    import numpy as np

    return np.array([profile_data.get(k, 0) for k in profile_channel_order()], dtype=float)


def run_cohort_analysis(
    input_file: Path,
    output_dir: Path,
    genome: str = "hg19",
    signatures_set: str = "COSMICv3",
    config: dict[str, bool] | None = None,
) -> dict[str, Any]:
    """Run comprehensive cohort analysis.

    Args:
        input_file: Path to input MAF/VCF file
        output_dir: Directory to store output files
        genome: Genome assembly (hg19, hg38, mm10)
        signatures_set: Signature set to use (COSMICv3, COSMICv2, KUCAB)
        config: Dict with analysis options (classify, cluster, hotspots, motifs)

    Returns:
        Dict with analysis results and paths to output files
    """
    if config is None:
        config = {}
    elif isinstance(config, str):
        # Handle case where config might be a JSON string
        try:
            config = json.loads(config)
        except (json.JSONDecodeError, TypeError):
            logger.warning(f"Invalid config format: {config}, using empty dict")
            config = {}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get genome file path from genome manager
    genome_manager = GenomeManager()
    genome_path = genome_manager.get_genome_path(genome)

    if not genome_path.exists():
        raise FileNotFoundError(
            f"Genome {genome} not found at {genome_path}. "
            f"Download it with: mutagene fetch genome {genome}"
        )

    results: dict[str, Any] = {
        "samples": 0,
        "mutations": 0,
        "profiles": {},
        "signatures": {},
        "classification": {},
        "files": [],
    }

    # Variants the caller rejected are excluded unless the run asks for them.
    keep_filtered = bool(config.get("keep_filtered", False))
    filter_column = config.get("filter_column") or "FILTER"

    # Extract tar.gz once; plain files pass through unchanged
    input_file = extract_input_file(Path(input_file), output_dir)

    try:
        # Step 1: Generate mutational profiles
        logger.info(f"Reading mutations from {input_file}")
        profile_file = output_dir / "profile.txt"

        with open_input_file(input_file, "rt") as infile:
            with open(profile_file, "w") as outfile:
                profile_stats = calc_profile(
                    [infile],
                    outfile,
                    str(genome_path),
                    fmt="auto",
                    keep_filtered=keep_filtered,
                    filter_column=filter_column,
                )

        # The mismatch count comes from the parser that does the comparison
        # rather than from counting its log records. The old handler watched a
        # logger around calc_profile, but the message it looked for is emitted
        # by a different code path and only after the handler was removed, so
        # it could never fire and wrong-assembly runs completed silently (#99).
        genome_warning = assembly_mismatch_warning(
            profile_stats, profile_stats.get("loaded", 0), genome
        )
        if genome_warning is not None:
            results["genome_warning"] = genome_warning

        # Parse the profile data for visualization
        profile_data = {}
        with open(profile_file) as pf:
            for line in pf:
                if line.strip():
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        mutation_type, count = parts
                        profile_data[mutation_type] = int(float(count))

        results["profiles"] = profile_data
        results["files"].append(
            {"type": "profile", "path": str(profile_file), "size": profile_file.stat().st_size}
        )

        # Step 2: Count mutations and samples
        # Use read_mutations to get per-sample data (needed for multi-sample support)
        from mutagene.io.context_window import read_mutations

        input_format = detect_input_format(input_file)
        logger.info(f"Reading per-sample mutations as {input_format}")

        with open_input_file(input_file, "rt") as infile:
            samples_mutations, _, processing_stats = read_mutations(
                input_format,
                infile,
                str(genome_path),
                window_size=1,
                keep_filtered=keep_filtered,
                filter_column=filter_column,
            )
            results["mutations"] = processing_stats.get("loaded", 0)
            # Detect number of samples from mutations dict
            num_samples = len(samples_mutations.keys()) if samples_mutations else 1
            results["samples"] = num_samples
            logger.info(f"Detected {num_samples} sample(s) in input file")

        # Step 3: Signature decomposition
        logger.info("Performing signature decomposition")
        try:
            import numpy as np

            from mutagene.io.profile import read_signatures
            from mutagene.profiles.profile import get_multisample_mutational_profile
            from mutagene.signatures.identify import decompose_mutational_profile_counts

            # Load signatures (use the signature set specified in config)
            W, signature_names = read_signatures(signatures_set, only=None)

            logger.info(f"Loaded {len(signature_names)} signatures from {signatures_set}")
            logger.info(f"Signature matrix shape: {W.shape}")

            # Convert profile_data dict to numpy array in signature-matrix row order.
            # NB: this must match W's row order, not alphabetical key order.
            profile_array = profile_dict_to_array(profile_data)

            logger.info(f"Profile array shape: {profile_array.shape}, sum: {profile_array.sum()}")

            # Decompose the pooled profile
            _, _, decomp_results = decompose_mutational_profile_counts(
                profile_array,
                (W, signature_names),
                func="MLE",
                others_threshold=0.0,
                enable_dummy=True,
            )

            logger.info(f"Decomposition returned {len(decomp_results)} results")

            # Format results for display
            signature_summary: dict[str, Any] = {
                "method": "MLE",
                "signatures": {},
                "total_mutations": int(sum(profile_data.values())),
                "per_sample": {},
            }

            # Extract signatures with non-zero contributions
            for sig_result in decomp_results:
                mutations = sig_result.get("mutations", 0)
                if (
                    mutations
                    and isinstance(mutations, (int, float, np.integer, np.floating))
                    and mutations > 0
                ):
                    sig_name = sig_result["name"]
                    exposure = sig_result["score"]
                    signature_summary["signatures"][sig_name] = {
                        "exposure": float(exposure),
                        "mutations": int(mutations),
                    }

            logger.info(
                f"Found {len(signature_summary['signatures'])} non-zero signatures in pooled profile"
            )

            # If multi-sample, decompose each sample individually
            if num_samples > 1:
                logger.info(f"Performing per-sample decomposition for {num_samples} samples")
                samples_profiles = get_multisample_mutational_profile(
                    samples_mutations, counts=True
                )

                for sample_id, sample_profile in samples_profiles.items():
                    # sample_profile is already a list of 96 floats in the correct order
                    sample_array = np.array(sample_profile, dtype=float)

                    # Decompose this sample
                    _, _, sample_decomp = decompose_mutational_profile_counts(
                        sample_array,
                        (W, signature_names),
                        func="MLE",
                        others_threshold=0.0,
                        enable_dummy=True,
                    )

                    # Extract non-zero signatures for this sample
                    sample_sigs = {}
                    for sig_result in sample_decomp:
                        mutations = sig_result.get("mutations", 0)
                        if (
                            mutations
                            and isinstance(mutations, (int, float, np.integer, np.floating))
                            and mutations > 0
                        ):
                            sig_name = sig_result["name"]
                            sample_sigs[sig_name] = {
                                "exposure": float(sig_result["score"]),
                                "mutations": int(mutations),
                            }

                    signature_summary["per_sample"][sample_id] = {
                        "signatures": sample_sigs,
                        "total_mutations": int(sample_array.sum()),
                    }

                logger.info(
                    f"Completed per-sample decomposition for {len(signature_summary['per_sample'])} samples"
                )

            results["signatures"] = signature_summary

            decomposition_file = output_dir / "decomposition.json"
            with open(decomposition_file, "w") as df:
                json.dump(signature_summary, df, indent=2)

            results["files"].append(
                {
                    "type": "decomposition",
                    "path": str(decomposition_file),
                    "size": decomposition_file.stat().st_size,
                }
            )
        except Exception as e:
            import traceback

            error_trace = traceback.format_exc()
            logger.error(f"Signature decomposition failed: {e}")
            logger.error(error_trace)
            results["signatures"] = {
                "method": "MLE",
                "signatures": {},
                "total_mutations": 0,
                "per_sample": {},
                "error": str(e),
            }

        # Step 4: Classification (if requested)
        if config.get("classify", False):
            logger.info("Running cancer type classification")
            try:
                classification_file = output_dir / "classification.json"
                with open(classification_file, "w") as cf:
                    json.dump({"predicted_type": "Unknown", "confidence": 0.0}, cf)

                results["files"].append(
                    {
                        "type": "classification",
                        "path": str(classification_file),
                        "size": classification_file.stat().st_size,
                    }
                )
            except Exception as e:
                logger.warning(f"Classification failed: {e}")

        # Step 5: Clustering (if requested and multi-sample)
        if config.get("cluster", False) and num_samples > 1:
            logger.info("Performing hierarchical clustering on signature exposures")
            try:
                import numpy as np
                from scipy.cluster.hierarchy import dendrogram, linkage
                from scipy.spatial.distance import pdist
                from sklearn.decomposition import PCA

                # Get per-sample signature data
                if results.get("signatures", {}).get("per_sample"):
                    per_sample_data = results["signatures"]["per_sample"]

                    # Build signature exposure matrix (samples x signatures)
                    signature_set = set()
                    for sample_data in per_sample_data.values():
                        signature_set.update(sample_data["signatures"].keys())

                    all_signatures = sorted(signature_set)
                    sample_ids = sorted(per_sample_data.keys())

                    # Create exposure matrix
                    exposure_rows = []
                    for sample_id in sample_ids:
                        sample_exposures = []
                        for sig_name in all_signatures:
                            exposure = (
                                per_sample_data[sample_id]["signatures"]
                                .get(sig_name, {})
                                .get("exposure", 0.0)
                            )
                            sample_exposures.append(exposure)
                        exposure_rows.append(sample_exposures)

                    exposure_matrix = np.array(exposure_rows)
                    logger.info(f"Exposure matrix shape: {exposure_matrix.shape}")

                    # Perform hierarchical clustering
                    if len(sample_ids) > MAX_CLUSTERING_SAMPLES:
                        logger.warning(
                            f"Skipping clustering: {len(sample_ids)} samples exceeds the "
                            f"limit of {MAX_CLUSTERING_SAMPLES} (pairwise distance matrix "
                            f"grows quadratically)"
                        )
                    elif len(sample_ids) >= 2:
                        # Calculate cosine distance
                        distances = pdist(exposure_matrix, metric="cosine")
                        linkage_matrix = linkage(distances, method="average")

                        # Generate dendrogram data
                        dend = dendrogram(linkage_matrix, no_plot=True, labels=sample_ids)

                        # PCA for 2D visualization
                        if exposure_matrix.shape[0] >= 2:
                            pca = PCA(
                                n_components=min(
                                    2, exposure_matrix.shape[0], exposure_matrix.shape[1]
                                )
                            )
                            pca_coords = pca.fit_transform(exposure_matrix)

                            # Prepare clustering results
                            clustering_data = {
                                "method": "hierarchical_average",
                                "distance_metric": "cosine",
                                "samples": sample_ids,
                                "dendrogram": {
                                    "icoord": dend["icoord"],
                                    "dcoord": dend["dcoord"],
                                    "leaves": dend["leaves"],
                                    "ivl": dend["ivl"],  # sample labels in dendrogram order
                                },
                                "pca_coordinates": {
                                    sample_id: {
                                        "x": float(pca_coords[i, 0]),
                                        "y": (
                                            float(pca_coords[i, 1])
                                            if pca_coords.shape[1] > 1
                                            else 0.0
                                        ),
                                    }
                                    for i, sample_id in enumerate(sample_ids)
                                },
                                "explained_variance": pca.explained_variance_ratio_.tolist(),
                            }

                            results["clustering"] = clustering_data

                            cluster_file = output_dir / "clustering.json"
                            with open(cluster_file, "w") as cf:
                                json.dump(clustering_data, cf, indent=2)

                            results["files"].append(
                                {
                                    "type": "clustering",
                                    "path": str(cluster_file),
                                    "size": cluster_file.stat().st_size,
                                }
                            )

                            logger.info(f"Clustering completed for {len(sample_ids)} samples")
                        else:
                            logger.warning("Not enough samples for PCA visualization")
                    else:
                        logger.warning("Need at least 2 samples for clustering")
                else:
                    logger.warning("No per-sample signature data available for clustering")

            except Exception as e:
                import traceback

                error_trace = traceback.format_exc()
                logger.error(f"Clustering failed: {e}")
                logger.error(error_trace)
        elif config.get("cluster", False) and num_samples == 1:
            logger.info("Skipping clustering: only 1 sample detected")

        # Step 6: Driver mutation ranking (hotspot detection)
        if config.get("hotspots", False):
            logger.info("Ranking driver mutations")
            try:
                import pandas as pd

                from mutagene.io.cohorts import read_cohort_mutations_from_tar
                from mutagene.io.protein_mutations_MAF import read_protein_mutations_MAF
                from mutagene.mutability.mutability import (
                    THRESHOLD_DRIVER,
                    THRESHOLD_PASSENGER,
                    rank,
                )

                # Read protein mutations
                with open_input_file(input_file, "rt") as infile:
                    protein_mutations, protein_stats = read_protein_mutations_MAF(
                        infile,
                        str(genome_path),
                        keep_filtered=keep_filtered,
                        filter_column=filter_column,
                    )

                logger.info(f"Loaded {protein_stats.get('loaded', 0)} protein mutations")

                if protein_mutations:
                    # Get profile for mutability calculation (signature-matrix row order)
                    profile_array_for_rank = profile_dict_to_array(profile_data)

                    # Try to load precalculated cohort data
                    cohorts_file = Path.home() / ".mutagene" / "cohorts.tar.gz"
                    cohort_profile = profile_array_for_rank
                    cohort_size = len(samples_mutations.keys()) if samples_mutations else 1
                    cohort_aa_mutations = None

                    # Check if cohorts file exists and try to use it
                    if cohorts_file.exists():
                        # Try common cancer types
                        cohort_name = config.get("cohort", None)
                        if cohort_name:
                            try:
                                cohort_profile, cohort_size, cohort_aa_mutations, _ = (
                                    read_cohort_mutations_from_tar(str(cohorts_file), cohort_name)
                                )
                                logger.info(
                                    f"Using precalculated cohort: {cohort_name} (N={cohort_size})"
                                )
                            except Exception as e:
                                logger.warning(f"Could not load cohort {cohort_name}: {e}")

                    # Run ranking and write to TSV file
                    drivers_file = output_dir / "drivers.tsv"
                    with open(drivers_file, "w") as df:
                        rank(
                            protein_mutations,
                            df,
                            cohort_profile,
                            cohort_aa_mutations,
                            cohort_size,
                            THRESHOLD_DRIVER,
                            THRESHOLD_PASSENGER,
                            provenance={
                                "mutagene_version": __version__,
                                "command": "rank",
                                "input_file": Path(input_file).name,
                                "genome": genome,
                                "profile_source": (
                                    f"precalculated cohort {config['cohort']}"
                                    if cohort_aa_mutations is not None
                                    else "input sample"
                                ),
                                "cohort_size": cohort_size,
                                "observed_mutations_source": (
                                    f"precalculated cohort {config['cohort']}"
                                    if cohort_aa_mutations is not None
                                    else "input sample"
                                ),
                                "threshold_driver": THRESHOLD_DRIVER,
                                "threshold_passenger": THRESHOLD_PASSENGER,
                            },
                        )

                    # Read the TSV back as JSON for storage. The provenance
                    # header is comment-prefixed, so it must be skipped here.
                    if drivers_file.exists() and drivers_file.stat().st_size > 0:
                        drivers_df = pd.read_csv(drivers_file, sep="\t", comment="#")
                        drivers_data = drivers_df.to_dict("records")

                        # Load known driver genes for annotation
                        known_drivers_file = (
                            Path(__file__).parent.parent / "data" / "known_drivers.json"
                        )
                        known_drivers = {}
                        if known_drivers_file.exists():
                            with open(known_drivers_file) as kdf:
                                known_drivers_json = json.load(kdf)
                                known_drivers = known_drivers_json.get("genes", {})
                            logger.info(f"Loaded {len(known_drivers)} known driver genes")

                        # Annotate with known driver information
                        for driver in drivers_data:
                            gene = driver.get("gene", "")
                            if gene in known_drivers:
                                driver["known_driver"] = True
                                driver["role"] = known_drivers[gene].get("role")
                                driver["tier"] = known_drivers[gene].get("tier")
                                driver["associated_cancers"] = known_drivers[gene].get(
                                    "cancers", []
                                )
                            else:
                                driver["known_driver"] = False

                        # Store in database
                        results["drivers"] = {
                            "total": len(drivers_data),
                            "drivers": [d for d in drivers_data if d.get("label") == "Driver"],
                            "potential_drivers": [
                                d for d in drivers_data if d.get("label") == "Potential driver"
                            ],
                            "passengers": [
                                d for d in drivers_data if d.get("label") == "Passenger"
                            ],
                            "known_driver_hits": [d for d in drivers_data if d.get("known_driver")],
                        }

                        logger.info(
                            f"Identified {len(results['drivers']['drivers'])} driver mutations"
                        )
                        logger.info(
                            f"Found {len(results['drivers']['known_driver_hits'])} mutations in known driver genes"
                        )

                        results["files"].append(
                            {
                                "type": "drivers",
                                "path": str(drivers_file),
                                "size": drivers_file.stat().st_size,
                            }
                        )
                else:
                    logger.warning("No protein mutations found for driver ranking")

            except Exception as e:
                import traceback

                error_trace = traceback.format_exc()
                logger.error(f"Driver ranking failed: {e}")
                logger.error(error_trace)

        # Step 7: Motif enrichment (if requested)
        if config.get("motifs", False):
            logger.info("Analyzing motif enrichment")
            try:
                motifs_file = output_dir / "motifs.json"
                with open(motifs_file, "w") as mf:
                    json.dump({"motifs": []}, mf)

                results["files"].append(
                    {"type": "motifs", "path": str(motifs_file), "size": motifs_file.stat().st_size}
                )
            except Exception as e:
                logger.warning(f"Motif enrichment failed: {e}")

        logger.info("Analysis completed successfully")
        return results

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise
