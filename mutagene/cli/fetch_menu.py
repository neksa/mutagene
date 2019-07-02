import sys
import logging
from mutagene.io.fetch import fetch_genome, fetch_cohorts, fetch_examples


logger = logging.getLogger(__name__)
genome_error_message = 'requires genome name argument -g hg19, hg38, mm10, see http://hgdownload.cse.ucsc.edu/downloads.html for more'


class FetchMenu(object):
    def __init__(self, parser):
        parser.description = 'Download data from remote repositories and API'
        # parser.add_argument('resource', choices=['cohorts', 'genome', 'GDC', 'MSKC', 'ICGC'])
        subparsers = parser.add_subparsers(
            dest='resource',
            title='subcommands',
            description='Choose data source',
            help='additional help available for subcommands')
        subparsers.required = True

        subparsers.add_parser('examples', description="""
This command will download sample1.maf and sample2.vcf files from the MutaGene server.
These files are only needed to reproduce examples in the manual.""")

        cohorts_parser = subparsers.add_parser('cohorts', help='cohorts help', description="""
This command will download mutational profiles and counts of observed mutations for cohorts available in online repositories.
The files are downloaded in one bundle cohorts.tar.gz that will be saved in the current directory.
Currently, only cohorts representing cancer types in COSMIC are provided.
Cohorts are required for ranking of mutations, because ranking relies upon counts of observed mutations and cancer type-specific profiles.
""")
        cohorts_parser_action = cohorts_parser.add_mutually_exclusive_group(required=True)
        cohorts_parser_action.add_argument('--list', '-l', action='store_true', help='List available cohorts')
        cohorts_parser_action.add_argument('--cohort', '-c', type=str, help='Specify cohort')
        cohorts_parser.add_argument('source', choices=('COSMIC', 'GDC', 'MSKC', 'ICGC'), default='COSMIC', help='Choose source data repository')

        genome_parser = subparsers.add_parser('genome', description="""
This command will download reference genome assembly sequence in 2bit format from the UCSC genome browser website.

You need to specify the name of genome assembly in --genome (-g) argument.

Partial download is supported: if the process is interupted run the same command again to continue downloading.
""")
        genome_parser.add_argument('--genome', '-g', type=str, help='hg38, hg19, mm10 according to UCSC genome browser nomenclature', required=True)

    @classmethod
    def examples(cls, args):
        # print('Cohorts', args)
        if args.resource == 'examples':
            fetch_examples()
            logger.info("Example files saved to current directory")

    @classmethod
    def cohorts(cls, args):
        # print('Cohorts', args)
        getattr(cls, 'cohorts_' + args.source)(args)

    @classmethod
    def cohorts_COSMIC(cls, args):
        fetch_cohorts()
        logger.info("cohorts.tar.gz saved to current directory")

    @classmethod
    def cohorts_GDC(cls, args):
        # fetch_cohorts()
        logger.warning("GDC currently not supported")

    @classmethod
    def cohorts_MSKC(cls, args):
        if args.list:
            import pandas as pd
            with pd.option_context('display.max_rows', None):  # , 'display.max_columns', None
                print(pd.read_csv("https://www.cbioportal.org/webservice.do?cmd=getCancerStudies", sep="\t"))
        elif args.cohort:
            logger.info("Will download cohort " + args.cohort)
        else:
            logger.warning("Cohort not specified, use --list to display all available data sets")

    @classmethod
    def cohorts_ICGC(cls, args):
        # fetch_cohorts()
        logger.warning("ICGC currently not supported")

    @classmethod
    def genome(cls, args):
        # print('Genome', args)
        if args.resource == 'genome':
            if not args.genome:
                logger.warning(genome_error_message)
                return
            fetch_genome(args.genome)
            logger.info("Twobit file saved to current directory")

    @classmethod
    def callback(cls, args):
        if not args.resource:
            logger.warning('No resource specified')
            sys.exit(1)

        # print('FetchMenu', args.resource)
        getattr(cls, args.resource)(args)
