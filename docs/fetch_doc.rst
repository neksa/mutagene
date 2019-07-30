==================================
Description
==================================
Use mutagene fetch to download example data files, genome assemblies, and cohorts to run mutagene commands

==================================
Fetch Subpackage Documentation
==================================

mutagene fetch {examples, cohorts, genome}

Fetch example mutational files (automatically download to current directory)

    Description: This command will download sample1.maf and sample2.vcf files from the MutaGene
    server.

    Command: run "mutagene fetch examples"

    Arguments:
     -h, --help  show this help message and exit

    Usage in MutaGene: These files are only needed to reproduce examples in the manual. Files can be inputted as
    infiles in other MutaGene subpackages.

    Data Fetched:
    Both sample1.maf and sample2.vcf are downloaded


Fetch cohorts

    Description: This command will download mutational profiles and counts of observed
    mutations for cohorts available in online repositories. The files are
    downloaded in one bundle cohorts.tar.gz that will be saved in the current
    directory. Currently, only cohorts representing cancer types in COSMIC are
    provided.


    Command: "mutagene fetch (--list | --cohort COHORT) {COSMIC,GDC,MSKC,ICGC}"

    Arguments:
     -h, --help  show this help message and exit
     -l, --list  List available cohorts
     -c, --cohort  Specify cohort

    Usage in MutaGene: Cohorts are required for ranking of mutations, because ranking
    relies upon counts of observed mutations and cancer type-specific profiles.


Fetch genome assemblies

    Description: This command will download reference genome assembly sequence in 2bit format
    from the UCSC genome browser website. You need to specify the name of genome
    assembly in --genome (-g) argument. Partial download is supported: if the
    process is interrupted run the same command again to continue downloading. Genome files are needed for mutagene motif,
    rank, profile, and signature commands.

    Command: "mutagene fetch genome -g <name of genome assembly>"

    Arguments:
     -h, --help         show this help message and exit
     -g, --genome
                        hg38, hg19, mm10 according to UCSC genome browser
                        nomenclature

    Usage in MutaGene:
    Genome assemblies required for genome argument in MutaGene subpackages.

    Examples:

    To download hg19 file:
      $ mutagene fetch genome -g hg19

    To download hg38 file:
     $ mutagene fetch genome -g hg38

    Data Fetched:
    Any reference genome available for download on the `UCSC Sequence and Annotation Data Webpage <http://hgdownload.soe.ucsc.edu/downloads.html>`_



