# Mutagene
Python library for mutational analysis


mutagene -i TCRBOA1-T-WEX_TCRBOA1-N-WEX.vcf -o out --signatures 5 --genome /Users/agoncear/data/hg38.2bit identify
mutagene identify -i TCRBOA1-T-WEX_TCRBOA1-N-WEX.vcf -o out --signatures 5 --genome /Users/agoncear/data/hg38.2bit

mutagene calc_profile --infile TCGA-50-6593-01.maf.txt --outfile TCGA-50-6593-01.profile --genome /Users/agoncear/data/hg38.2bit

mutagene rank -m aa -g /Users/agoncear/data/hg38.2bit -i TCGA-50-6593-01.maf.txt
