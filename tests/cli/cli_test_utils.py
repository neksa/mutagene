import hashlib, sys
from mutagene.__main__ import MutaGeneApp


TEST_DIR = 'test-reports'
COHORTS_FILE = 'cohorts.tar.gz'
TEST_FILE_MAP = { COHORTS_FILE: 'https://www.ncbi.nlm.nih.gov/research/mutagene/static/data/cohorts.tar.gz',
                    'sample1.maf': 'https://www.ncbi.nlm.nih.gov/research/mutagene/static/data/sample1.maf',
                    'hg19.2bit': 'https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit' }
                    # 'chrY.fa.gz': 'https://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrY.fa.gz'


# Method to run MutaGene CLI commands in a manner similar to actual CLI execution
def run_with_args(cmd, cmd_args):
    # Copy original sys.argv from pytest execution
    argv_orig = sys.argv.copy()

    del sys.argv[1:]

    sys.argv.append('-v')
    sys.argv.append(cmd)

    for a in cmd_args:
        sys.argv.append(a)

    # Execute MutaGene with current arguments
    MutaGeneApp()

    # Restore pytest sys.argv
    sys.argv = argv_orig


# Method to return a number of lines of a file
def get_file_lines(outfile, num_lines):
    out_lines = []
    out_fh = open(outfile, 'r')

    for i in range(num_lines):
        out_lines.append(out_fh.readline())
        if i >= num_lines: break

    out_fh.close()

    return out_lines


# https://stackoverflow.com/a/3431838
def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()
