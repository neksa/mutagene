"""Genome reference management for webapp."""

import logging
import threading
from pathlib import Path

from mutagene.io.fetch import download_from_url

logger = logging.getLogger(__name__)


class GenomeManager:
    """Manage genome reference files for the webapp."""

    SUPPORTED_GENOMES = ["hg19", "hg38", "mm10", "mm9"]
    GENOME_URLS = {
        "hg19": "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit",
        "hg38": "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit",
        "mm10": "https://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit",
        "mm9": "https://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit",
    }

    def __init__(self, genomes_dir: Path | None = None):
        """Initialize genome manager.

        Args:
            genomes_dir: Directory to store genome files.
                         Defaults to ~/.mutagene/genomes/
        """
        if genomes_dir is None:
            genomes_dir = Path.home() / ".mutagene" / "genomes"

        self.genomes_dir = Path(genomes_dir)
        self.genomes_dir.mkdir(parents=True, exist_ok=True)

        # Serializes downloads per assembly: two concurrent requests for the
        # same genome would otherwise write the same partial file.
        self._download_locks: dict[str, threading.Lock] = {}
        self._locks_guard = threading.Lock()

    def _lock_for(self, genome: str) -> threading.Lock:
        with self._locks_guard:
            return self._download_locks.setdefault(genome, threading.Lock())

    def get_genome_path(self, genome: str) -> Path:
        """Get the path to a genome file.

        Args:
            genome: Genome assembly name (e.g., 'hg19', 'hg38')

        Returns:
            Path to the .2bit file

        Raises:
            ValueError: If *genome* is not a supported assembly name. Rejecting
                unknown names here keeps request-supplied values from being
                interpolated into a filesystem path.
        """
        if genome not in self.SUPPORTED_GENOMES:
            raise ValueError(
                f"Unsupported genome: {genome!r}. Supported: {', '.join(self.SUPPORTED_GENOMES)}"
            )
        return self.genomes_dir / f"{genome}.2bit"

    def is_downloaded(self, genome: str) -> bool:
        """Check if a genome has been downloaded.

        Args:
            genome: Genome assembly name

        Returns:
            True if genome file exists. Unknown assembly names return False
            rather than raising, so callers can probe freely.
        """
        try:
            return self.get_genome_path(genome).exists()
        except ValueError:
            return False

    def get_available_genomes(self) -> list[str]:
        """Get list of available (downloaded) genomes.

        Returns:
            List of genome names that are available
        """
        return [g for g in self.SUPPORTED_GENOMES if self.is_downloaded(g)]

    def get_missing_genomes(self) -> list[str]:
        """Get list of supported but not downloaded genomes.

        Returns:
            List of genome names that need to be downloaded
        """
        return [g for g in self.SUPPORTED_GENOMES if not self.is_downloaded(g)]

    def download_genome(self, genome: str, progress_callback=None) -> bool:
        """Download a genome reference file.

        Args:
            genome: Genome assembly name
            progress_callback: Optional callback function for progress updates

        Returns:
            True if download successful, False otherwise
        """
        if genome not in self.GENOME_URLS:
            logger.error(f"Unsupported genome: {genome}")
            return False

        with self._lock_for(genome):
            return self._download_locked(genome)

    def _download_locked(self, genome: str) -> bool:
        # Another thread may have finished this download while we waited.
        if self.is_downloaded(genome):
            return True

        url = self.GENOME_URLS[genome]
        final_path = self.get_genome_path(genome)
        # Download to a sibling temp file and rename only on success, so an
        # interrupted download never leaves a truncated file that
        # is_downloaded() would report as a usable genome.
        partial_path = final_path.with_name(final_path.name + ".part")

        try:
            logger.info(f"Downloading {genome} from {url}")
            download_from_url(url, str(partial_path))
            partial_path.replace(final_path)
            logger.info(f"Successfully downloaded {genome} to {final_path}")
            return True
        except Exception as e:
            logger.error(f"Failed to download {genome}: {e}")
            # Clean up partial download
            partial_path.unlink(missing_ok=True)
            return False

    def check_and_download_required_genomes(
        self, required: list[str] | None = None, auto_download: bool = False
    ) -> dict:
        """Check for required genomes and optionally download them.

        Args:
            required: List of required genomes. Defaults to ['hg19', 'hg38']
            auto_download: If True, automatically download missing genomes

        Returns:
            Dict with 'available', 'missing', and 'downloaded' keys
        """
        if required is None:
            required = ["hg19", "hg38"]

        available = []
        missing = []
        downloaded = []

        for genome in required:
            if self.is_downloaded(genome):
                available.append(genome)
            else:
                missing.append(genome)
                if auto_download:
                    logger.info(f"Auto-downloading missing genome: {genome}")
                    if self.download_genome(genome):
                        downloaded.append(genome)
                        available.append(genome)

        return {"available": available, "missing": missing, "downloaded": downloaded}

    def get_genome_info(self) -> dict:
        """Get information about all supported genomes.

        Returns:
            Dict mapping genome names to info dicts
        """
        info = {}
        for genome in self.SUPPORTED_GENOMES:
            path = self.get_genome_path(genome)
            info[genome] = {
                "downloaded": path.exists(),
                "path": str(path),
                "size_mb": round(path.stat().st_size / 1024 / 1024, 1) if path.exists() else None,
                "url": self.GENOME_URLS.get(genome),
            }
        return info
