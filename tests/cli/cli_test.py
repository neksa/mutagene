import os, shutil

from mutagene.__main__ import MutaGeneApp
from tests.cli import cli_test_utils


def test_fetch():
    # mutagene -v fetch cohorts MSKCC --cohort paac_jhu_2014
    cli_test_utils.run_with_args('fetch', ['cohorts', 'MSKCC', '--cohort', 'paac_jhu_2014'])

    file_name = 'paac_jhu_2014.tar.gz'
    os.rename(f'./{file_name}', f'{cli_test_utils.TEST_DIR}/{file_name}')

    file_md5sum = cli_test_utils.md5sum(f'{cli_test_utils.TEST_DIR}/{file_name}')

    assert file_md5sum == 'b7709f55eaeade1b1c6102d134b16c18'


def test_motif(artifactory_circleci):
    infile = f'{cli_test_utils.TEST_DIR}/sample1.maf'
    outfile = f'{cli_test_utils.TEST_DIR}/cli-motif-sample1.txt'
    genome = f'{cli_test_utils.TEST_DIR}/hg19.2bit'

    # mutagene -v motif -i sample1.maf -g hg19 --motif "C[A>T]" --strand A -o test-reports/cli-motif-sample1.txt
    cli_test_utils.run_with_args('motif', ['-i', infile, '-o', outfile, '-g', genome, '--motif', 'C[A>T]', '--strand', 'A'])

    out_lines = []
    in_fh = open(outfile, 'r')

    for i in range(2):
        out_lines.append(in_fh.readline())
        if i >= 2: break

    in_fh.close()

    assert out_lines[1].startswith('TCGA-50-6593-01A-11D-1753-08\tCustom motif\tC[A>T]\tany strand')


def test_profile(artifactory_circleci):
    infile = 'tests/motifs/data/vcf/data.vcf'
    outfile = f'{cli_test_utils.TEST_DIR}/cli-profile-test.txt'
    genome = f'{cli_test_utils.TEST_DIR}/hg19.2bit'  # 'tests/motifs/data/test_genome.2bit'

    # mutagene -v profile -i sample1.maf -g hg19 -o test-reports/cli-profile-sample1.txt
    cli_test_utils.run_with_args('profile', ['-i', infile, '-o', outfile, '-g', genome])

    out_lines = []
    in_fh = open(outfile, 'r')

    for i in range(3):
        out_lines.append(in_fh.readline())
        if i >= 3: break

    in_fh.close()

    assert out_lines[0] == 'A[C>A]A\t0\n'
    assert out_lines[1] == 'A[C>G]A\t2\n'
    assert out_lines[2] == 'A[C>T]A\t4\n'


def test_rank(artifactory_circleci):
    infile = f'{cli_test_utils.TEST_DIR}/sample1.maf'
    outfile = f'{cli_test_utils.TEST_DIR}/cli-rank-sample1-pancancer.txt'
    genome = f'{cli_test_utils.TEST_DIR}/hg19.2bit'

    cp_cohorts = False
    if not os.path.isfile(f'./{cli_test_utils.COHORTS_FILE}'):
        cp_cohorts = True
        shutil.copyfile(f'{cli_test_utils.TEST_DIR}/{cli_test_utils.COHORTS_FILE}', f'./{cli_test_utils.COHORTS_FILE}')

    # mutagene -v rank -g hg19 -i sample1.maf -c pancancer -o test-reports/cli-rank-sample1-pancancer.txt
    cli_test_utils.run_with_args('rank', ['-i', infile, '-o', outfile, '-g', genome, '-c', 'pancancer'])

    if cp_cohorts == True:
        os.remove(f'./{cli_test_utils.COHORTS_FILE}')

    out_lines = []
    in_fh = open(outfile, 'r')

    for i in range(2):
        out_lines.append(in_fh.readline())
        if i >= 2: break

    in_fh.close()

    assert out_lines[1].startswith('CPXM2\tT536M\t')


def test_signature(artifactory_circleci):
    infile = f'{cli_test_utils.TEST_DIR}/sample1.maf'
    outfile = f'{cli_test_utils.TEST_DIR}/cli-signature-sample1.txt'
    genome = f'{cli_test_utils.TEST_DIR}/hg19.2bit'

    # mutagene -v signature identify -i sample1.maf -g hg19 -s5 -o test-reports/cli-signature-sample1.txt
    cli_test_utils.run_with_args('signature', ['identify', '-i', infile, '-o', outfile, '-g', genome, '-s5'])

    out_lines = []
    in_fh = open(outfile, 'r')

    for i in range(2):
        out_lines.append(in_fh.readline())
        if i >= 2: break

    in_fh.close()

    assert out_lines[1].startswith('TCGA-50-6593-01A-11D-1753-08\t1\t0.5361\t182')
