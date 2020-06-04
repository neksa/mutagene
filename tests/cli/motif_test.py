from tests.cli import cli_test_utils


DEFAULT_INFILE = f'{cli_test_utils.TEST_DIR}/sample1.maf'
DEFAULT_GENOME = f'{cli_test_utils.TEST_DIR}/hg19.2bit'


# mutagene -v motif -i sample1.maf -g hg19 --motif "C[A>T]" --strand A -o test-reports/cli-motif-sample1-CA_T.txt
def test_motif(artifactory_circleci):
    infile = DEFAULT_INFILE
    outfile = f'{cli_test_utils.TEST_DIR}/cli-motif-sample1-CA_T.txt'
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args('motif', ['-i', infile, '-o', outfile, '-g', genome, '--motif', 'C[A>T]', '--strand', 'A'])

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    assert out_lines[1].startswith('TCGA-50-6593-01A-11D-1753-08\tCustom motif\tC[A>T]\tany strand')


# mutagene -v motif -i sample1.maf -g hg19 --strand A -o test-reports/cli-motif-sample1.txt
def test_motif_5_1(artifactory_circleci):
    infile = DEFAULT_INFILE
    outfile = f'{cli_test_utils.TEST_DIR}/cli-motif-sample1.txt'
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args('motif', ['-i', infile, '-o', outfile, '-g', genome, '--strand', 'A'])

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    assert out_lines[1].startswith('TCGA-50-6593-01A-11D-1753-08\tC>T in CpG\t[C>T]G\tany strand\t')


# mutagene -v motif -i sample1.maf -g hg19 --strand A -w 20 -o test-reports/cli-motif-sample1-w20.txt
def test_motif_5_3(artifactory_circleci):
    infile = DEFAULT_INFILE
    outfile = f'{cli_test_utils.TEST_DIR}/cli-motif-sample1-w20.txt'
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args('motif', ['-i', infile, '-o', outfile, '-g', genome, '--strand', 'A', '-w', '20'])

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    assert out_lines[1].startswith('TCGA-50-6593-01A-11D-1753-08\tAPOBEC3G\tC[C>K]R\tany strand\t')
    assert out_lines[2].startswith('TCGA-50-6593-01A-11D-1753-08\tC>T in CpG\t[C>T]G\tany strand\t')
