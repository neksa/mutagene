import json
import os
import pandas as pd

import logging
logger = logging.getLogger(__name__)


def write_motif_matches(outfile, motif_matches):
    """Convert motif matches list of dictionaries to CSV format"""
    if len(motif_matches) > 0:
        header = motif_matches[0].keys()
        d = pd.DataFrame.from_records(motif_matches, columns=header).sort_values(by='pvalue')
        d.to_csv(outfile, sep="\t", index=False, float_format='%g')


def get_known_motifs():
    """Load a pre-bundled set of motifs from data directory"""
    dirname = os.path.dirname(os.path.realpath(__file__))
    fname = dirname + "/../data/motifs/motifs.json"
    with open(fname, 'r') as f:
        motifs = json.load(f)
    return motifs
