"""Genome reference management for webapp."""

import logging
from pathlib import Path
from typing import List, Optional

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

    def __init__(self, genomes_dir: Optional[Path] = None):
        """Initialize genome manager.

        Args:
            genomes_dir: Directory to store genome files.
                         Defaults to ~/.mutagene/genomes/
        """
        if genomes_dir is None:
            genomes_dir = Path.home() / ".mutagene" / "genomes"

        self.genomes_dir = Path(genomes_dir)
        self.genomes_dir.mkdir(parents=True, exist_ok=True)

    def get_genome_path(self, genome: str) -> Path:
        """Get the path to a genome file.

        Args:
            genome: Genome assembly name (e.g., 'hg19', 'hg38')

        Returns:
            Path to the .2bit file
        """
        return self.genomes_dir / f"{genome}.2bit"

    def is_downloaded(self, genome: str) -> bool:
        """Check if a genome has been downloaded.

        Args:
            genome: Genome assembly name

        Returns:
            True if genome file exists
        """
        return self.get_genome_path(genome).exists()

    def get_available_genomes(self) -> List[str]:
        """Get list of available (downloaded) genomes.

        Returns:
            List of genome names that are available
        """
        return [g for g in self.SUPPORTED_GENOMES if self.is_downloaded(g)]

    def get_missing_genomes(self) -> List[str]:
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

        url = self.GENOME_URLS[genome]
        dst = str(self.get_genome_path(genome))

        try:
            logger.info(f"Downloading {genome} from {url}")
            download_from_url(url, dst)
            logger.info(f"Successfully downloaded {genome} to {dst}")
            return True
        except Exception as e:
            logger.error(f"Failed to download {genome}: {e}")
            # Clean up partial download
            if Path(dst).exists():
                Path(dst).unlink()
            return False

    def check_and_download_required_genomes(
        self, required: Optional[List[str]] = None, auto_download: bool = False
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
