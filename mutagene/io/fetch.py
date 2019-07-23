import os
# import io
import sys
import requests

from tqdm import tqdm

import logging
logger = logging.getLogger(__name__)


# ########### URL ##############
def download_from_url(url, dst):
    """
    Download file with progressbar

    @param: url to download file
    @param: dst place to put the file
    """
    chunk_size = 1024  # 1048576

    try:
        file_size = int(requests.head(url).headers["Content-Length"])
    except:
        logger.warning("Could not access URL " + url)
        sys.exit(1)

    first_byte = os.path.getsize(dst) if os.path.exists(dst) else 0

    if first_byte >= file_size:
        logger.warning("Looks like the file has been downloaded already: " + dst)
        return file_size

    header = {"Range": "bytes=%s-%s" % (first_byte, file_size)}
    try:
        r = requests.get(url, headers=header, stream=True)
        file_size = int(r.headers['Content-Length'])
    except:
        logger.warning("Could not access URL " + url)
        sys.exit(1)

    def generate(raw, chunk_size):
        while True:
            chunk = raw.read(chunk_size)
            if not chunk:
                break
            yield chunk

    # need to access raw stream because r.iter_content() deflates .gz automatically
    # opening file in append mode
    with open(dst, 'ab') as fp:
        for chunk in tqdm(
                generate(r.raw, chunk_size=chunk_size),
                initial=first_byte // chunk_size,
                total=file_size // chunk_size,
                unit='KB',
                desc=dst,
                leave=True):
            fp.write(chunk)
    return file_size


def fetch_genome(name):
    # ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
    # http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
    assembly = name.rsplit('.', 1)[0] if '.2bit' in name else name
    assembly = assembly.lower()
    mode = 'http'

    url = f'{mode}://hgdownload.cse.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}.2bit'
    dst = f'{assembly}.2bit'
    download_from_url(url, dst)


def fetch_cohorts():
    url = "https://www.ncbi.nlm.nih.gov/research/mutagene/static/data/cohorts.tar.gz"
    dst = "cohorts.tar.gz"
    download_from_url(url, dst)


def fetch_examples():
    url = "https://www.ncbi.nlm.nih.gov/research/mutagene/static/data/"
    files = ["sample1.maf", "sample2.vcf"]
    for fname in files:
        download_from_url(url + fname, fname)


def fetch_MSKCC(dataset):
    """ Download dataset from cBioPortal https://www.cbioportal.org/datasets """
    url = "http://download.cbioportal.org/{}.tar.gz".format(dataset)
    dst = "{}.tar.gz".format(dataset)
    download_from_url(url, dst)
