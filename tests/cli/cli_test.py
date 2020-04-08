import unittest
import hashlib, os, requests, shutil, sys
from mutagene.__main__ import MutaGeneApp


TEST_DIR = 'test-reports'
COHORTS_FILE = 'cohorts.tar.gz'
TEST_FILE_LIST = [COHORTS_FILE, 'sample1.maf', 'hg19.2bit']


def setup_module(module):
    ARTIFACTORY_ROOT_URL = 'http://169.53.172.72:8081/artifactory/generic-local/mutagene'
    ARTIFACTORY_USER = 'mutagene'
    ARTIFACTORY_PASSWD = 'w8$X2:Eb[Ug7Di6@'

    for f in TEST_FILE_LIST:
        if not os.path.isfile(f'{TEST_DIR}/{f}'):
            if os.path.isfile(f'./{f}'):
                shutil.copyfile(f'./{f}', f'{TEST_DIR}/{f}')
            else:
                r = requests.get(f'{ARTIFACTORY_ROOT_URL}/{f}', auth=(ARTIFACTORY_USER, ARTIFACTORY_PASSWD))
                outfile = open(f'{TEST_DIR}/{f}', 'wb')
                outfile.write(r.content)
                outfile.close()


def teardown_module(module):
    if 'CIRCLECI' in os.environ and os.environ['CIRCLECI'] == 'true':
        for f in TEST_FILE_LIST:
            if os.path.isfile(f'{TEST_DIR}/{f}'):
                os.remove(f'{TEST_DIR}/{f}')


def set_argv(cmd, cmd_args):
    argv_orig = sys.argv

    del sys.argv[1:]

    sys.argv.append('-v')
    sys.argv.append(cmd)

    for a in cmd_args:
        sys.argv.append(a)

    return argv_orig


# https://stackoverflow.com/a/3431838
def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            hash.update(block)
    return hash.hexdigest()


class CliTestCases(unittest.TestCase):
    def test_fetch(self):
        # mutagene -v fetch cohorts MSKCC --cohort paac_jhu_2014
        argv_orig = set_argv('fetch', ['cohorts', 'MSKCC', '--cohort', 'paac_jhu_2014'])
        MutaGeneApp()
        sys.argv = argv_orig

        file_name = 'paac_jhu_2014.tar.gz'
        os.rename(f'./{file_name}', f'{TEST_DIR}/{file_name}')

        file_md5sum = md5sum(f'{TEST_DIR}/{file_name}')

        assert file_md5sum == 'b7709f55eaeade1b1c6102d134b16c18'


    def test_motif(self):
        infile = f'{TEST_DIR}/sample1.maf'
        outfile = f'{TEST_DIR}/cli-motif-sample1.txt'
        genome = f'{TEST_DIR}/hg19.2bit'

        # mutagene -v motif -i sample1.maf -g hg19 --motif "C[A>T]" --strand A -o test-reports/cli-motif-sample1.txt
        argv_orig = set_argv('motif', ['-i', infile, '-o', outfile, '-g', genome, '--motif', 'C[A>T]', '--strand', 'A'])
        MutaGeneApp()
        sys.argv = argv_orig

        out_lines = []
        in_fh = open(outfile, 'r')

        for i in range(2):
            out_lines.append(in_fh.readline())
            if i >= 2: break

        in_fh.close()

        assert out_lines[1].startswith('TCGA-50-6593-01A-11D-1753-08\tCustom motif\tC[A>T]\tany strand')


    def test_profile(self):
        infile = 'tests/motifs/data/vcf/data.vcf'
        outfile = f'{TEST_DIR}/cli-profile-test.txt'
        genome = f'{TEST_DIR}/hg19.2bit'  # 'tests/motifs/data/test_genome.2bit'

        # mutagene -v profile -i sample1.maf -g hg19 -o test-reports/cli-profile-sample1.txt
        argv_orig = set_argv('profile', ['-i', infile, '-o', outfile, '-g', genome])
        MutaGeneApp()
        sys.argv = argv_orig

        out_lines = []
        in_fh = open(outfile, 'r')

        for i in range(3):
            out_lines.append(in_fh.readline())
            if i >= 3: break

        in_fh.close()

        assert out_lines[0] == 'A[C>A]A\t0\n'
        assert out_lines[1] == 'A[C>G]A\t2\n'
        assert out_lines[2] == 'A[C>T]A\t4\n'


    def test_rank(self):
        infile = f'{TEST_DIR}/sample1.maf'
        outfile = f'{TEST_DIR}/cli-rank-sample1-pancancer.txt'
        genome = f'{TEST_DIR}/hg19.2bit'

        cp_cohorts = False
        if not os.path.isfile(f'./{COHORTS_FILE}'):
            cp_cohorts = True
            shutil.copyfile(f'{TEST_DIR}/{COHORTS_FILE}', f'./{COHORTS_FILE}')

        # mutagene -v rank -g hg19 -i sample1.maf -c pancancer -o test-reports/cli-rank-sample1-pancancer.txt
        argv_orig = set_argv('rank', ['-i', infile, '-o', outfile, '-g', genome, '-c', 'pancancer'])
        MutaGeneApp()
        sys.argv = argv_orig

        if cp_cohorts == True:
            os.remove(f'./{COHORTS_FILE}')

        out_lines = []
        in_fh = open(outfile, 'r')

        for i in range(2):
            out_lines.append(in_fh.readline())
            if i >= 2: break

        in_fh.close()

        assert out_lines[1].startswith('CPXM2\tT536M\t')


    def test_signature(self):
        infile = f'{TEST_DIR}/sample1.maf'
        outfile = f'{TEST_DIR}/cli-signature-sample1.txt'
        genome = f'{TEST_DIR}/hg19.2bit'

        # mutagene -v signature identify -i sample1.maf -g hg19 -s5 -o test-reports/cli-signature-sample1.txt
        argv_orig = set_argv('signature', ['identify', '-i', infile, '-o', outfile, '-g', genome, '-s5'])
        MutaGeneApp()
        sys.argv = argv_orig

        out_lines = []
        in_fh = open(outfile, 'r')

        for i in range(2):
            out_lines.append(in_fh.readline())
            if i >= 2: break

        in_fh.close()

        assert out_lines[1].startswith('TCGA-50-6593-01A-11D-1753-08\t1\t0.5361\t182')


if __name__ == '__main__':
#    print(os.getcwd() + '\n' + sys.argv)
    unittest.main()
#    print(sys.argv)
