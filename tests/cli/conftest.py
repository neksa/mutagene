import pytest
import os, requests, shutil

from tests.cli import cli_test_utils


# Download test files from public sources if not present
@pytest.fixture(scope='module')
def artifactory_circleci():
    # Ensure TEST_DIR exists
    os.makedirs(cli_test_utils.TEST_DIR, exist_ok=True)

    # Download test files from public sources
    for filename, url in cli_test_utils.TEST_FILE_MAP.items():
        target_path = f'{cli_test_utils.TEST_DIR}/{filename}'

        if not os.path.isfile(target_path):
            # File may already exist at execution root, e.g., when run by CI
            # If so, copy file to test-reports for this set of tests
            if os.path.isfile(f'./{filename}'):
                shutil.copyfile(f'./{filename}', target_path)
            # Otherwise, download from public source
            else:
                print(f"Downloading {filename} from {url}")
                # Use stream=True to avoid automatic decompression
                # Important for .tar.gz files that need to stay compressed
                r = requests.get(url, stream=True)
                r.raise_for_status()
                r.raw.decode_content = False  # Keep content as-is, don't decode gzip
                with open(target_path, 'wb') as outfile:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            outfile.write(chunk)

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
