"""Analysis pipeline integration for webapp."""

import gzip
import json
import logging
from pathlib import Path
from typing import Any, Dict

from mutagene.profiles.profile import calc_profile

from .genome_manager import GenomeManager

logger = logging.getLogger(__name__)


def extract_input_file(file_path: Path, output_dir: Path) -> Path:
    """Extract tar.gz archives to a stable location; return path to the mutation file.

    For non-archive files the original path is returned unchanged.
    Extracted files go into *output_dir* so they are cleaned up with the analysis.
    """
    import tarfile

    if str(file_path).endswith(".tar.gz") or str(file_path).endswith(".tgz"):
        with tarfile.open(file_path, "r:gz") as tar:
            maf_candidates = [
                m
                for m in tar.getmembers()
                if m.isfile() and ("mutation" in m.name.lower() or m.name.endswith(".maf"))
            ]
            if not maf_candidates:
                raise ValueError("No mutation file found in tar.gz archive")

            maf_member = maf_candidates[0]
            logger.info(f"Extracting {maf_member.name} from tarball")

            extract_dir = output_dir / "extracted"
            extract_dir.mkdir(exist_ok=True)
            extracted_path = (extract_dir / Path(maf_member.name).name).resolve()
            if not str(extracted_path).startswith(str(extract_dir.resolve())):
                raise ValueError(
                    f"Tar member {maf_member.name} would extract outside target directory"
                )
            tar.extract(maf_member, path=extract_dir)
            return extracted_path

    return file_path


def open_input_file(file_path: Path, mode: str = "rt"):
    """Open a mutation file, handling gzip transparently."""
    if str(file_path).endswith(".gz"):
        return gzip.open(file_path, mode, encoding="utf-8")
    return open(file_path, mode, encoding="utf-8")


def run_cohort_analysis(
    input_file: Path,
    output_dir: Path,
    genome: str = "hg19",
    signatures_set: str = "COSMICv3",
    config: Dict[str, bool] = None,
) -> Dict[str, Any]:
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

    results = {
        "samples": 0,
        "mutations": 0,
        "profiles": {},
        "signatures": {},
        "classification": {},
        "files": [],
    }

    # Extract tar.gz once; plain files pass through unchanged
    input_file = extract_input_file(Path(input_file), output_dir)

    # Track genome mismatch warnings
    genome_mismatch_count = 0

    try:
        # Step 1: Generate mutational profiles
        logger.info(f"Reading mutations from {input_file}")
        profile_file = output_dir / "profile.txt"

        # Count genome mismatches by intercepting log messages from context_window.
        # Coupled to the warning string in io/context_window.py — update if that changes.
        import logging as _logging

        class MismatchCounter(_logging.Handler):
            def __init__(self):
                super().__init__()
                self.mismatch_count = 0

            def emit(self, record):
                if "REF allele does not match" in record.getMessage():
                    self.mismatch_count += 1

        mismatch_handler = MismatchCounter()
        ctx_logger = _logging.getLogger("mutagene.io.context_window")
        ctx_logger.addHandler(mismatch_handler)

        try:
            with open_input_file(input_file, "rt") as infile:
                with open(profile_file, "w") as outfile:
                    calc_profile([infile], outfile, str(genome_path), fmt="auto")
        finally:
            ctx_logger.removeHandler(mismatch_handler)

        genome_mismatch_count = mismatch_handler.mismatch_count

        if genome_mismatch_count > 0:
            logger.warning(f"Found {genome_mismatch_count} genome assembly mismatches")
            if genome_mismatch_count > 20:
                logger.error(
                    f"High number of genome mismatches ({genome_mismatch_count}) suggests wrong assembly!"
                )
                results["genome_warning"] = {
                    "mismatch_count": genome_mismatch_count,
                    "message": f"Found {genome_mismatch_count} reference allele mismatches. This suggests the wrong genome assembly was selected. Please try {('hg38' if genome == 'hg19' else 'hg19')}.",
                }

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

        with open_input_file(input_file, "rt") as infile:
            samples_mutations, _, processing_stats = read_mutations(
                "MAF", infile, str(genome_path), window_size=1
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

            # Convert profile_data dict to numpy array in correct order
            profile_keys = sorted(profile_data.keys())
            profile_array = np.array([profile_data[k] for k in profile_keys], dtype=float)

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
            signature_summary = {
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
                    all_signatures = set()
                    for sample_data in per_sample_data.values():
                        all_signatures.update(sample_data["signatures"].keys())

                    all_signatures = sorted(list(all_signatures))
                    sample_ids = sorted(per_sample_data.keys())

                    # Create exposure matrix
                    exposure_matrix = []
                    for sample_id in sample_ids:
                        sample_exposures = []
                        for sig_name in all_signatures:
                            exposure = (
                                per_sample_data[sample_id]["signatures"]
                                .get(sig_name, {})
                                .get("exposure", 0.0)
                            )
                            sample_exposures.append(exposure)
                        exposure_matrix.append(sample_exposures)

                    exposure_matrix = np.array(exposure_matrix)
                    logger.info(f"Exposure matrix shape: {exposure_matrix.shape}")

                    # Perform hierarchical clustering
                    if len(sample_ids) >= 2:
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
                        infile, str(genome_path)
                    )

                logger.info(f"Loaded {protein_stats.get('loaded', 0)} protein mutations")

                if protein_mutations:
                    # Get profile for mutability calculation
                    profile_array_for_rank = np.array(
                        [profile_data[k] for k in sorted(profile_data.keys())], dtype=float
                    )

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
                        )

                    # Read the TSV back as JSON for storage
                    if drivers_file.exists() and drivers_file.stat().st_size > 0:
                        drivers_df = pd.read_csv(drivers_file, sep="\t")
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
