
* pytest 
 - parametrizing tests
 - fixtures

* meaningful readme.md with
 - main section
 - Installation
 - Changes
 - Documentation
 - Usage (printout of argparse usage)
 
* readthedocs documentation

* clustergrammer
* rpy2

* cell lines examples

- AL LUAD data analysis
- Nitin LUAD analysis
- benchmark clean up
- threshold?
- packaging & CI & PIP
- BayesianOptimization



mutagene rank -g hg38.2bit -i sample1.maf -o ranking.txt


mutagene -i TCRBOA1-T-WEX_TCRBOA1-N-WEX.vcf -o out --signatures 5 --genome /Users/agoncear/data/hg38.2bit identify
mutagene identify -i TCRBOA1-T-WEX_TCRBOA1-N-WEX.vcf -o out --signatures 5 --genome /Users/agoncear/data/hg38.2bit

mutagene calc_profile --infile TCGA-50-6593-01.maf.txt --outfile TCGA-50-6593-01.profile --genome /Users/agoncear/data/hg38.2bit

mutagene rank -m aa -g /Users/agoncear/data/hg38.2bit -i TCGA-50-6593-01.maf.txt



twobitreader.download.save_genome(name, destdir=None, mode='ftp')[source]
tries to download a genome from UCSC by name

for example, ‘hg19’ is at ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit



https://plastid.readthedocs.io/en/latest/quickstart.html

https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/

https://github.com/mskcc/vcf2maf
https://github.com/PoisonAlien/maftools
https://riptutorial.com/bioinformatics/example/14312/mutation-annotation-format--maf-

http://docs.h5py.org/en/latest/high/dataset.html#dataset