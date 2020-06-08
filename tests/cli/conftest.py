import pytest
import os, requests, shutil

from tests.cli import cli_test_utils


# Download large test articles from Artifactory if not present in current environment
@pytest.fixture(scope='module')
def artifactory_circleci(store):
    # Start of module setup code, using pytest-securestore
    ARTIFACTORY_ROOT_URL = store['artifactory']['url']
    ARTIFACTORY_USER = store['artifactory']['username']
    ARTIFACTORY_PASSWD = store['artifactory']['password']

    for f in cli_test_utils.TEST_FILE_MAP.keys():
        if not os.path.isfile(f'{cli_test_utils.TEST_DIR}/{f}'):
            # File may already exist at execution root, e.g., when run by CircleCI
            # If so, copy file to test-reports for this set of tests
            if os.path.isfile(f'./{f}'):
                shutil.copyfile(f'./{f}', f'{cli_test_utils.TEST_DIR}/{f}')
            # Otherwise, first attempt to download from Artifactory
            else:
                r = requests.get(f'{ARTIFACTORY_ROOT_URL}/{f}', auth=(ARTIFACTORY_USER, ARTIFACTORY_PASSWD))
                # If request failed, download from original source specified in map
                try:
                    r.raise_for_status()
                    outfile = open(f'{cli_test_utils.TEST_DIR}/{f}', 'wb')
                    outfile.write(r.content)
                    outfile.close()
                except:
                    r = requests.get(cli_test_utils.TEST_FILE_MAP[f])
                    outfile = open(f'{cli_test_utils.TEST_DIR}/{f}', 'wb')
                    outfile.write(r.content)
                    outfile.close()

                    # Attempt to upload this file to Artifactory
                    # https://github.com/jfrog/project-examples/blob/master/bash-example/deploy-file.sh
                    # file_md5sum = md5sum(f'{cli_test_utils.TEST_DIR}/{f}')
                    r = requests.put(f'{ARTIFACTORY_ROOT_URL}/{f}', auth=(ARTIFACTORY_USER, ARTIFACTORY_PASSWD),
                            files = {'file': open(f'{cli_test_utils.TEST_DIR}/{f}', 'rb')})  # headers = { 'X-Checksum-Md5': file_md5sum }

    # Copy COHORTS_FILE to ./ for use in mutagene rank tests
    cp_cohorts = False
    if not os.path.isfile(f'./{cli_test_utils.COHORTS_FILE}'):
        cp_cohorts = True
        shutil.copyfile(f'{cli_test_utils.TEST_DIR}/{cli_test_utils.COHORTS_FILE}', f'./{cli_test_utils.COHORTS_FILE}')

    # End of module setup code
    yield
    # Start of teardown code

    # If COHORTS_FILE was copied to ./, remove it
    if cp_cohorts == True:
        os.remove(f'./{cli_test_utils.COHORTS_FILE}')

    # Teardown code cleans up test articles from TEST_DIR if running in CircleCI to prevent saving as artifacts
    if 'CIRCLECI' in os.environ and os.environ['CIRCLECI'] == 'true':
        for f in cli_test_utils.TEST_FILE_MAP.keys():
            if os.path.isfile(f'{cli_test_utils.TEST_DIR}/{f}'):
                os.remove(f'{cli_test_utils.TEST_DIR}/{f}')
