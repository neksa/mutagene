"""Resolution of the ``--genome`` (``-g``) argument shared by the analysis subcommands.

``mutagene fetch genome`` downloads assemblies into ``~/.mutagene/genomes``, so a
bare assembly name is resolved there before falling back to the current directory.

The genomes directory lives here rather than in ``mutagene.webapp.genome_manager``
because importing anything from ``mutagene.webapp`` pulls in Flask, an optional
dependency the CLI has to work without.
"""

import os
from pathlib import Path

DEFAULT_GENOMES_DIR = Path.home() / ".mutagene" / "genomes"


class GenomeNotFoundError(FileNotFoundError):
    """Raised when a genome assembly cannot be located in any known location."""


def _twobit_name(genome: str) -> str:
    return genome if genome.endswith(".2bit") else genome + ".2bit"


def resolve_genome(genome: str, genomes_dir: Path | None = None) -> str:
    """Resolve a ``-g`` value to a usable genome assembly location.

    Args:
        genome: Either a path to a 2bit file or a bare assembly name (e.g. ``hg19``).
        genomes_dir: Directory holding downloaded assemblies, defaults to
            :data:`DEFAULT_GENOMES_DIR`.

    Returns:
        The value to hand to the readers: an existing path when one was found,
        otherwise the original value unchanged (preserving the historical
        current-directory behaviour).

    Raises:
        GenomeNotFoundError: If no assembly was found in any searched location.
    """
    if genomes_dir is None:
        genomes_dir = DEFAULT_GENOMES_DIR

    # An existing file path is used as given
    if os.path.isfile(genome):
        return genome

    twobit = _twobit_name(genome)

    # A bare assembly name resolves against the directory fetch downloads into
    if os.sep not in genome and (os.altsep is None or os.altsep not in genome):
        downloaded = genomes_dir / twobit
        if downloaded.is_file():
            return str(downloaded)

    # Historical behaviour: <name>.2bit relative to the current directory
    if os.path.isfile(twobit):
        return genome

    raise GenomeNotFoundError(
        f"Genome assembly '{genome}' was not found. Searched:\n"
        f"    {genomes_dir / twobit}\n"
        f"    {Path.cwd() / twobit}\n"
        f"Download it with: mutagene fetch genome -g {genome.rsplit('.2bit', 1)[0]}"
    )
