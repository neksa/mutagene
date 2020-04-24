from tests.cli import cli_test_utils


DEFAULT_INFILE = f'{cli_test_utils.TEST_DIR}/sample1.maf'
DEFAULT_GENOME = f'{cli_test_utils.TEST_DIR}/hg19.2bit'


# mutagene -v rank -g hg19 -i sample1.maf -c pancancer -o test-reports/cli-rank-sample1-pancancer.txt
def test_rank(artifactory_circleci):
    infile = DEFAULT_INFILE
    outfile = f'{cli_test_utils.TEST_DIR}/cli-rank-sample1-pancancer.txt'
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args('rank', ['-i', infile, '-o', outfile, '-g', genome, '-c', 'pancancer'])

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    assert out_lines[1].startswith('CPXM2\tT536M\t')


# mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas -o test-reports/cli-rank-sample1-gcb_lymphomas.txt
def test_rank_4_1(artifactory_circleci):
    infile = DEFAULT_INFILE
    outfile = f'{cli_test_utils.TEST_DIR}/cli-rank-sample1-gcb_lymphomas.txt'
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args('rank', ['-i', infile, '-o', outfile, '-g', genome, '-c', 'gcb_lymphomas'])

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    assert out_lines[1].startswith('BOD1L\tT2810S\t')


# mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas -tp 0.0003 -o test-reports/cli-rank-sample1-gcb_lymphomas-tp3.txt
def test_rank_4_2(artifactory_circleci):
    infile = DEFAULT_INFILE
    outfile = f'{cli_test_utils.TEST_DIR}/cli-rank-sample1-gcb_lymphomas-tp3.txt'
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args('rank', ['-i', infile, '-o', outfile, '-g', genome, '-c', 'gcb_lymphomas', '-tp', '0.0003'])

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    assert out_lines[1].startswith('WNT8B\tR231C\t')


# mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas -td 0.000009 -o test-reports/cli-rank-sample1-gcb_lymphomas-td9.txt
def test_rank_4_3(artifactory_circleci):
    infile = DEFAULT_INFILE
    outfile = f'{cli_test_utils.TEST_DIR}/cli-rank-sample1-gcb_lymphomas-td9.txt'
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args('rank', ['-i', infile, '-o', outfile, '-g', genome, '-c', 'gcb_lymphomas', '-td', '0.000009'])

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    assert out_lines[1].startswith('C1orf69\tE244V\t')


# mutagene rank -i sample1.maf -g hg19 -c gcb_lymphomas -n 20 -o test-reports/cli-rank-sample1-gcb_lymphomas-n20.txt
def test_rank_4_4(artifactory_circleci):
    infile = DEFAULT_INFILE
    outfile = f'{cli_test_utils.TEST_DIR}/cli-rank-sample1-gcb_lymphomas-n20.txt'
    genome = DEFAULT_GENOME

    cli_test_utils.run_with_args('rank', ['-i', infile, '-o', outfile, '-g', genome, '-c', 'gcb_lymphomas', '-n', '20'])

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    assert out_lines[1].startswith('BOD1L\tT2810S\t')
