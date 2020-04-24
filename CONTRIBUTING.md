# Contributing

## Setup for General Development

* Fork and clone this repo, e.g. using the instructions at [GitHub: Fork a repo](https://help.github.com/en/github/getting-started-with-github/fork-a-repo)
  * Navigate to [neksa/mutagene](https://github.com/neksa/mutagene), click the Fork button at the top right
  * `git clone https://github.com/<your_github_username>/mutagene.git`
  * `cd mutagene`
  * `git remote add upstream https://github.com/neksa/mutagene.git`
  * Keep your fork in sync using the instructions at [GitHub: Syncing a fork](https://help.github.com/en/articles/syncing-a-fork)
* Use your package manager of choice to [install pipenv](https://pipenv.readthedocs.io/en/stable/install/#installing-pipenv)
* Create a new [pipenv](https://pipenv.readthedocs.io/en/stable/) environment
  * e.g., if you need to specify a `python` executable:  `pipenv install --python=/usr/bin/python3.7`
* Install the package in development mode using [setup.py](https://setuptools.readthedocs.io/en/latest/)
  * `pipenv run python setup.py develop`

## Adding Features or Fixing Bugs

* Create and check out a feature or bugfix branch
  * e.g., `git checkout -b feature/motif_extension`
* Add and commit your changes
* Update `README.md` when needed
* Add a description of your changes to `CHANGELOG.md`
* Ensure all tests are passing, as described below
* [Submit a pull request](https://help.github.com/en/articles/about-pull-requests) to the upstream repo with a description of your changes
* Ensure the new branch is mergeable

## Testing

* Write additional tests to cover new features as necessary
* Please make sure that all tests pass when running [pytest](https://docs.pytest.org/en/stable/)
* Our tests make use of the [pytest-securestore](https://pypi.org/project/pytest-securestore/) plugin
  to store credentials for connecting to an instance of [JFrog Artifactory](https://jfrog.com/open-source/#artifactory)
  that stores test articles.  In order to run tests using `pytest`, it is necessary to set the following
  environment variables:
  * `SECURE_STORE_FILE=./pytest-secure.yml.aes`
  * `SECURE_STORE_PASSWORD=`*contact [Alexander Goncearenco](mailto:gonceare@ncbi.nlm.nih.gov) or
    [Steven Handy](mailto:steven.handy@queensu.ca) for instructions*
  * To add a test article to Artifactory, update the `TEST_FILE_MAP` dictionary in the file `tests/cli/cli_test_utils.py`
  * If you wish to use your own Artifactory repository, the test scripts can obtain the articles from their original
    sources and upload them to your Artifactory repository.  See `pytest-secure.yml.template` for the required format
    for a file to contain your connection information.
  * To encrypt your YAML file, use the following commands in a `pipenv run python` interpreter:
      * `import pyAesCrypt`
      * `pyAesCrypt.encryptFile('./<your_artifactory_creds_file>.yml', './<your_artifactory_creds_file>.yml.aes', '<your_encryption_password>', (64 * 1024))`
* Run the tests using `pytest`
  * `pipenv run pytest --secure-store-filename=$SECURE_STORE_FILE --secure-store-password=$SECURE_STORE_PASSWORD tests`
