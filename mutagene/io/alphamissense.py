"""Cache and query the canonical AlphaMissense variant tables."""

import hashlib
import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path

import requests

logger = logging.getLogger(__name__)
RECORD = "https://zenodo.org/records/10813168"
DEFAULT_CACHE = Path.home() / ".mutagene" / "alphamissense"
SOURCES = {
    "hg19": {"size": 622293310, "md5": "3f14fba08c60b09a90a2ae14b7802e4c"},
    "hg38": {"size": 642961469, "md5": "9fd167735f16a1b87da6eb3e4c25fcb5"},
}
COLUMNS = (
    "chrom",
    "pos",
    "ref",
    "alt",
    "genome",
    "uniprot_id",
    "transcript_id",
    "protein_variant",
    "am_pathogenicity",
    "am_class",
)


def require_pysam():
    """Load the optional tabix dependency."""
    try:
        import pysam
    except ImportError as error:
        raise ValueError(
            'Install tabix support with pip install "mutagene[annotations]"'
        ) from error
    return pysam


def file_hash(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch_alphamissense(genome, cache_dir=None):
    """Download, verify, and index the pinned bgzip table without loading it into memory."""
    import fcntl

    pysam = require_pysam()
    if genome not in SOURCES:
        raise ValueError("AlphaMissense supports hg19 and hg38")
    directory = Path(cache_dir or DEFAULT_CACHE) / genome
    directory.mkdir(parents=True, exist_ok=True)
    filename = f"AlphaMissense_{genome}.tsv.gz"
    destination = directory / filename
    expected = SOURCES[genome]
    url = f"https://storage.googleapis.com/dm_alphamissense/{filename}"
    with (directory / ".download.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if destination.exists():
            if (
                destination.stat().st_size != expected["size"]
                or file_hash(destination, "md5") != expected["md5"]
            ):
                raise ValueError(
                    f"Cached AlphaMissense checksum differs: {destination}. Remove the corrupt file and retry."
                )
        else:
            partial = destination.with_suffix(destination.suffix + ".part")
            offset = partial.stat().st_size if partial.exists() else 0
            if offset > expected["size"]:
                raise ValueError(f"Partial download exceeds the expected size: {partial}")
            if offset < expected["size"]:
                logger.warning(f"Downloading {filename} ({expected['size'] / 1e6:.1f} MB)")
                with requests.get(
                    url, headers={"Range": f"bytes={offset}-"}, stream=True, timeout=(15, 120)
                ) as response:
                    response.raise_for_status()
                    if response.status_code == 206:
                        if not response.headers.get("Content-Range", "").startswith(
                            f"bytes {offset}-"
                        ):
                            raise ValueError("Server returns an unexpected download range")
                        mode = "ab"
                    elif response.status_code == 200:
                        mode = "wb"
                    else:
                        raise ValueError(f"Unexpected download status: {response.status_code}")
                    response.raw.decode_content = False
                    with partial.open(mode) as output:
                        while block := response.raw.read(1024 * 1024):
                            output.write(block)
            if (
                partial.stat().st_size != expected["size"]
                or file_hash(partial, "md5") != expected["md5"]
            ):
                raise ValueError(f"AlphaMissense download fails checksum verification: {partial}")
            os.replace(partial, destination)
        index = Path(str(destination) + ".tbi")
        if not index.exists():
            logger.warning("Creating the AlphaMissense tabix index")
            pysam.tabix_index(
                str(destination), seq_col=0, start_col=1, end_col=1, zerobased=False, force=True
            )
        metadata = {
            "genome": genome,
            "source_record": RECORD,
            "download_url": url,
            "license": "CC-BY-4.0",
            "citation_doi": "10.1126/science.adg7492",
            "size": destination.stat().st_size,
            "md5": expected["md5"],
            "sha256": file_hash(destination),
            "index_sha256": file_hash(index),
            "verified_at": datetime.now(timezone.utc).isoformat(),
            "coordinate_system": "1-based",
            "transcripts": "canonical GENCODE V27 (hg19) or V32 (hg38)",
        }
        manifest = Path(str(destination) + ".provenance.json")
        temporary = manifest.with_suffix(".tmp")
        temporary.write_text(json.dumps(metadata, indent=2) + "\n")
        os.replace(temporary, manifest)
    return destination


class AlphaMissenseLookup:
    """Match genome, allele, transcript, and protein change before returning a score."""

    def __init__(self, path, genome):
        pysam = require_pysam()
        if genome not in SOURCES:
            raise ValueError("Specify hg19 or hg38 for AlphaMissense lookup")
        self.path = Path(path)
        self.genome = genome
        self.metadata = {
            "path": str(self.path.resolve()),
            "genome": genome,
            "source_record": RECORD,
        }
        manifest = Path(str(path) + ".provenance.json")
        if manifest.exists():
            self.metadata.update(json.loads(manifest.read_text()))
            if self.metadata["genome"] != genome:
                raise ValueError("AlphaMissense assembly differs from the requested assembly")
        self.table = pysam.TabixFile(str(path))

    def close(self):
        self.table.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def lookup(self, chrom, pos, ref, alt, transcript, protein_variant):
        """Return a score only for an unambiguous missense match."""
        missing = {
            "alphamissense_score": None,
            "alphamissense_class": None,
            "alphamissense_status": "missing",
        }
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
            return dict(missing, alphamissense_status="unsupported_allele")
        protein_variant = protein_variant.removeprefix("p.")
        if not protein_variant or protein_variant[-1] in "*=X":
            return dict(missing, alphamissense_status="not_missense")
        chrom = str(chrom)
        chrom = chrom if chrom.startswith("chr") else "chr" + chrom
        if chrom not in self.table.contigs:
            return missing
        records = []
        for line in self.table.fetch(chrom, int(pos) - 1, int(pos)):
            fields = line.split("\t")
            if len(fields) != len(COLUMNS):
                raise ValueError("AlphaMissense record has an invalid number of columns")
            row = dict(zip(COLUMNS, fields))
            if row.get("genome") != self.genome:
                raise ValueError("AlphaMissense record uses a different genome assembly")
            if (row.get("pos"), row.get("ref"), row.get("alt")) == (str(pos), ref, alt):
                records.append(row)
        if not records:
            return missing
        # A versioned transcript must match exactly. Unversioned matches remain explicit.
        exact = "." in transcript
        matched = [
            r
            for r in records
            if (
                r["transcript_id"] == transcript
                if exact
                else r["transcript_id"].split(".")[0] == transcript
            )
        ]
        if not matched:
            return dict(missing, alphamissense_status="transcript_mismatch")
        matched = [r for r in matched if r["protein_variant"] == protein_variant]
        if not matched:
            return dict(missing, alphamissense_status="protein_mismatch")
        if len(matched) != 1:
            return dict(missing, alphamissense_status="ambiguous")
        row = matched[0]
        score = float(row["am_pathogenicity"])
        if not 0 <= score <= 1:
            raise ValueError("AlphaMissense score lies outside [0, 1]")
        return {
            "alphamissense_score": score,
            "alphamissense_class": row["am_class"],
            "alphamissense_status": "matched" if exact else "matched_unversioned_transcript",
            "alphamissense_transcript": row["transcript_id"],
            "alphamissense_uniprot": row["uniprot_id"],
        }


def annotate_protein_mutations(mutations, lookup):
    """Attach separate functional evidence without changing B-scores or driver labels."""
    annotations = {}
    for (gene, change), properties in mutations.items():
        variants = sorted(properties.get("genomic_variants", ()))
        matches = [
            lookup.lookup(*variant, properties["transcript"], change) for variant in variants
        ]
        if not matches:
            annotation = {
                "alphamissense_score": None,
                "alphamissense_class": None,
                "alphamissense_status": "missing_coordinates",
            }
        elif all(match == matches[0] for match in matches):
            annotation = matches[0]
        else:
            annotation = {
                "alphamissense_score": None,
                "alphamissense_class": None,
                "alphamissense_status": "ambiguous_genomic_mapping",
            }
        annotations[(gene, change)] = annotation
    return annotations
