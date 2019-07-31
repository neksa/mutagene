import pandas as pd

import logging
logger = logging.getLogger(__name__)


def write_motif_matches(outfile, motif_matches):
    if len(motif_matches) > 0:
        header = motif_matches[0].keys()
        d = pd.DataFrame.from_records(motif_matches, columns=header)
        d.to_csv(outfile, sep="\t", index=False, float_format='%g')
