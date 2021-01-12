import os

from tests.cli import cli_test_utils


def test_fetch():
    # mutagene -v fetch cohorts MSKCC --cohort paac_jhu_2014
    cli_test_utils.run_with_args('fetch', ['cohorts', 'MSKCC', '--cohort', 'paac_jhu_2014'])

    file_name = 'paac_jhu_2014.tar.gz'
    os.rename(f'./{file_name}', f'{cli_test_utils.TEST_DIR}/{file_name}')

    # Removed due to inconsistent md5 sums across repeated downloads, affecting both this test and CircleCI
    #file_md5sum = cli_test_utils.md5sum(f'{cli_test_utils.TEST_DIR}/{file_name}')
    #assert file_md5sum == 'acbf8c569c2b8f5684ccfb1e036743f0'

    file_size = os.path.getsize(f'{cli_test_utils.TEST_DIR}/{file_name}')
    assert file_size == 271744


def test_profile(artifactory_circleci):
    infile = 'tests/motifs/data/vcf/data.vcf'
    outfile = f'{cli_test_utils.TEST_DIR}/cli-profile-test.txt'
    genome = f'{cli_test_utils.TEST_DIR}/hg19.2bit'  # 'tests/motifs/data/test_genome.2bit'

    # mutagene -v profile -i sample1.maf -g hg19 -o test-reports/cli-profile-sample1.txt
    cli_test_utils.run_with_args('profile', ['-i', infile, '-o', outfile, '-g', genome])

    out_lines = cli_test_utils.get_file_lines(outfile, 3)

    assert out_lines[0] == 'A[C>A]A\t0\n'
    assert out_lines[1] == 'A[C>G]A\t2\n'
    assert out_lines[2] == 'A[C>T]A\t4\n'


def test_signature(artifactory_circleci):
    infile = f'{cli_test_utils.TEST_DIR}/sample1.maf'
    outfile = f'{cli_test_utils.TEST_DIR}/cli-signature-sample1.txt'
    genome = f'{cli_test_utils.TEST_DIR}/hg19.2bit'

    # mutagene -v signature identify -i sample1.maf -g hg19 -s5 -o test-reports/cli-signature-sample1.txt
    cli_test_utils.run_with_args('signature', ['identify', '-i', infile, '-o', outfile, '-g', genome, '-s5'])

    out_lines = cli_test_utils.get_file_lines(outfile, 2)

    assert out_lines[1].startswith('TCGA-50-6593-01A-11D-1753-08\tMGA-1\t0.532422\t180')
