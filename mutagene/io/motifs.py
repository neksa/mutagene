import csv

import logging
logger = logging.getLogger(__name__)


def write_motif_matches(outfile, motif_matches):
    if len(motif_matches) > 0:
        header = motif_matches[0].keys()
        writer = csv.DictWriter(outfile, fieldnames=header, dialect='excel-tab')
        writer.writeheader()
        for data in motif_matches:
            writer.writerow(data)
