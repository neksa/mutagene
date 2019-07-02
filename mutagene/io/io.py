# import json
# import requests

# import twobitreader as tbr

from collections import defaultdict
# from collections import namedtuple
from tqdm import tqdm

from mutagene.dna import nucleotides, complementary_nucleotide
# from mutagene.dna import codon_table
from mutagene.dna import chromosome_name_mapping

import logging
logger = logging.getLogger(__name__)
