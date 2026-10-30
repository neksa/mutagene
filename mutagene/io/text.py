"""Open plain or gzip-compressed mutation input for CLI commands."""

import argparse
import gzip
import shutil
import sys
import tempfile


def open_text_input(path):
    """Return a text stream for a path or standard input."""
    if path == "-":
        # Ranking reads the input twice; spool large streams to disk.
        stream = tempfile.SpooledTemporaryFile(
            max_size=8 * 1024 * 1024, mode="w+t", encoding="utf-8"
        )
        shutil.copyfileobj(sys.stdin, stream)
        stream.seek(0)
        return stream
    try:
        if str(path).lower().endswith(".gz"):
            return gzip.open(path, "rt", encoding="utf-8")
        return open(path, encoding="utf-8")
    except OSError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def detect_mutation_format(stream):
    """Identify VCF metadata or a MAF header without consuming the input."""
    if not stream.seekable():
        raise ValueError("Use --input-format for non-seekable input")
    position = stream.tell()
    try:
        for line in stream:
            if line.startswith(("##fileformat=VCF", "#CHROM")):
                return "VCF"
            if line.strip() and not line.startswith("#"):
                return "MAF"
    finally:
        stream.seek(position)
    return "MAF"
