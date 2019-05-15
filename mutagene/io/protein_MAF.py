# import os
import csv
# import numpy as np
import twobitreader as tbr

# from collections import defaultdict
from collections import namedtuple
# from itertools import cycle
from tqdm import tqdm

from mutagene.dna import nucleotides, complementary_nucleotide
from mutagene.dna import codon_table
# from mutagene.dna import chromosome_name_mapping

import logging
logger = logging.getLogger(__name__)


def read_MAF_with_genomic_context_fname(fname, genome, motifs=False):
    with open(fname) as infile:
        return read_MAF_with_genomic_context(infile, genome)


def read_MAF_with_genomic_context(infile, genome, motifs=False):
    mutations = []
    processing_stats = {'loaded': 0, 'skipped': 0, 'format': 'Unknown'}

    if not infile:
        logger.warning("No input file")
        return mutations, processing_stats

    if genome.upper() != 'MAF':
        fname = genome if genome.endswith('.2bit') else genome + '.2bit'
        twobit = tbr.TwoBitFile(fname)

    try:
        reader = csv.reader((row for row in infile if not row.startswith('#')), delimiter='\t')
        # get names from column headers
        header = next(reader)
        # MAF = namedtuple("MAF", map(str.lower, next(reader)))
        header = tuple(map(lambda s: s.replace('.', '_'), header))
        # print(header)
        MAF = namedtuple("MAF", header)
    except ValueError:
        raise
        logger.warning("MAF format not recognized")
        return mutations, processing_stats

    N_loaded = N_skipped = 0

    samples = set()

    for data in tqdm(map(MAF._make, reader)):
        # print(data)
        try:
            if hasattr(data, 'Tumor_Sample_Barcode'):
                samples.add(data.Tumor_Sample_Barcode)

            HGVSc = None
            if hasattr(data, 'cDNA_Change'):
                HGVSc = data.cDNA_Change
            if hasattr(data, 'TxChange'):
                HGVSc = data.TxChange

            if HGVSc is None:
                logger.warning("Could not find cDNA_Change or TxChange fields in MAF file")
                break
            if not HGVSc.startswith("c."):
                N_skipped += 1
                logger.debug("HGVSc c.")
                continue

            if '_' in HGVSc or 'del' in HGVSc:
                N_skipped += 1
                logger.debug("Skip indel")
                continue

            if '>' in HGVSc:
                # c.368C>T
                cDNA_position = int(HGVSc.split('.')[1][:-3])
            else:
                # c.C232T
                cDNA_position = int(HGVSc.split('.')[1][1:-1])
            offset = ((cDNA_position % 3) - 1) % 3  # determine codon offset: 0, 1, 2

            HGVSp = None
            if hasattr(data, 'Protein_Change'):
                HGVSp = data.Protein_Change
            if hasattr(data, 'HGVSp_Short'):
                HGVSp = data.HGVSp_Short
            if hasattr(data, 'AAChange'):
                HGVSp = data.AAChange

            if HGVSp is None:
                logger.warning("Could not find HGVSp_Short or Protein_Change or AAChange fields in MAF file")
                break

            if not HGVSp.startswith("p."):
                N_skipped += 1
                logger.debug("HGVSp p.")
                continue

            if '_' in HGVSp or HGVSp.endswith('fs'):
                N_skipped += 1
                logger.debug("HGVSp fs")
                continue

            protein_mutation = HGVSp.split('.')[1].upper()
            P = protein_mutation[0]
            Q = protein_mutation[-1]

            if not hasattr(data, 'Hugo_Symbol'):
                logger.warning("Could not find Hugo_Symbol in MAF file")
                break

            if not hasattr(data, 'Variant_Classification'):
                logger.warning("Could not find Variant_Classification in MAF file")
                break

            if data.Variant_Classification.lower() not in ("nonsense_mutation", "missense_mutation", "silent"):
                N_skipped += 1
                logger.debug("Variant_Classification")
                continue

            chrom = None
            start = None

            if hasattr(data, 'Chromosome'):
                chrom = data.Chromosome
            if hasattr(data, 'Start_position'):
                start = int(data.Start_position)
            if hasattr(data, 'Start_Position'):
                start = int(data.Start_Position)

            if not chrom or not start:
                logger.warning("Could not find Chromosome and Start_Position or Start_position in MAF file")
                break

            # if not hasattr(data, 'Codon_Change') and not hasattr(data, 'Codons'):
            #     logger.warning("MutaGen requires Codon_Change or Codons fields in MAF files")
            #     break

            # if hasattr(data, 'Codon_Change'):
            #     c1, c2 = data.Codon_Change.upper().split(")")[1].split(">")

            # if hasattr(data, 'Codons'):
            #     c1, c2 = data.Codons.upper().split("/")

            # if c1 == c2:
            #     N_skipped += 1
            #     logger.debug("Codon1 == Codon2")
            #     continue
            # for offset in range(3):
            #     if c1[offset] != c2[offset]:
            #         break

            chromosome = chrom if chrom.startswith('chr') else 'chr' + chrom

            if not motifs and hasattr(data, 'ref_context'):
                # context bundled in MAF file
                context = data.ref_context.upper()
                start = 10 + 1  # ref_context: 10 + 1 + 10
                seq5_fwd = context[start - offset - 2: start - offset + 3]
                seq5_rev_mirror = context[start - (2 - offset) - 2: start - (2 - offset) + 3]
                seq5_rev = "".join([complementary_nucleotide[x] for x in reversed(seq5_rev_mirror)])
            else:
                if genome.upper() == 'MAF':
                    logger.warning("ref_context not found in MAF file. Provide genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more")
                    break

                # try:
                #     window_size = 50
                #     seq_window = twobit[chromosome][start - window_size: start + window_size + 1]
                #     assert len(seq_window) == window_size * 2 + 1
                #     seq_window = seq_window.upper()

                #     seq_with_coords = list(zip(
                #         cycle([chrom]),
                #         range(start - window_size - 1, start + window_size),
                #         seq_window,
                #         cycle("+")))

                #     # print(seq_with_coords)
                # except:
                #     pass

                try:
                    seq5_fwd = twobit[chromosome][start - offset - 2: start - offset + 3].upper()
                    seq5_rev_mirror = twobit[chromosome][start - (2 - offset) - 2: start - (2 - offset) + 3].upper()
                except ValueError:
                    # could not read reference genome with given coordinates
                    logger.debug("Detected mismatch of mutated reference allele with the reference genome, check if using correct genome: " + data.Protein_Change)
                    # logger.warning("Chr {} {} {}".format(chromosome, start, offset))
                    N_skipped += 1
                    continue

                if len(set(seq5_fwd + seq5_rev_mirror) - set(nucleotides)) > 0:
                    # detected non-nucleotide letters (including N)
                    N_skipped += 1
                    logger.debug("Detected non-nucleotide letters (including N)")
                    continue
                seq5_rev = "".join([complementary_nucleotide[x] for x in reversed(seq5_rev_mirror)])

            if codon_table[seq5_fwd[1:4]] == P:
                seq5 = seq5_fwd
            elif codon_table[seq5_rev[1:4]] == P:
                seq5 = seq5_rev
            else:
                # print(seq5_fwd[1:4], seq5_rev[1:4], c1)
                N_skipped += 1
                # logger.debug(data.Reference_Allele)
                # logger.debug(data.Tumor_Seq_Allele1)
                # logger.debug(data.Tumor_Seq_Allele2)
                # logger.debug(context[:10] + " " + context[10] + " " + context[11:] + " " + seq5_fwd + " " + seq5_rev + " " + c1 + " " + c2)
                logger.debug("Sequence missmatch of mutation with the genome, check if using the correct genome assembly: " + HGVSp)
                continue
            # data.Chromosome,
            # data.Start_position,
            # data.Strand,
            # data.Codon_Change,

            N_loaded += 1
            mutations.append({
                'gene': data.Hugo_Symbol,
                'mutation': HGVSp.split(".")[1],
                'seq5': seq5,
            })
        except Exception as e:
            N_skipped += 1
            logger.debug("General MAF parsing exception " + str(e))
            continue
        # print(data.Codon_Change, data.Strand, offset, seq5)

        # print(data.cDNA_position, data.CDS_position, data.Protein_position)
        # print(data.Codon_Change)  #, data.CDS_position, data.Protein_position)
        # cDNA_position   Relative position of base pair in the cDNA sequence as a fraction. A "-" symbol is displayed as the numerator if the variant does not appear in cDNA
        # 54 - CDS_position   Relative position of base pair in coding sequence. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence
        # 55 - Protein_position
        # print(data.HGVSp)
        # Hugo_Symbol
        # Chromosome
        # Start_Position
        # End_Position
        # Strand
        # HGVSp
    processing_stats = {'loaded': N_loaded, 'skipped': N_skipped, 'format': 'MAF'}
    return mutations, processing_stats
