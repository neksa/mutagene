import requests
import json
import csv
import glob

# MUTAGENE_URL = "https://www.ncbi.nlm.nih.gov/research/mutagene"
#MUTAGENE_URL = "https://dev.ncbi.nlm.nih.gov/research/mutagene"
MUTAGENE_URL = "http://mwebdev2/research/mutagene"
# MUTAGENE_URL = "http://localhost:5000"

path = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/mutagene/data/skin/mutations/FI_D1_D2_combined_WGS.maf.txt"

def get_motifs(fname, assembly=37):
    """
    Identify mutational motifs in a VCF or MAF sample
    """
    url = MUTAGENE_URL + '/pub/api/identify/motifs'
    files = {'file': open(fname, 'rb')}
    r = requests.post(url, files=files, data={'assembly': assembly})
    if r.status_code == 200:
        return r.json()['motifs']


def print_motifs(motifs):
    """
    Printing the results of motif identification
    """
    if motifs is None:
        print("Empty")
        return
    for m in motifs:
        print("{}\t{}\t{:.2f}\t{:.2e}\t{}".format(
            m['name'], m['motif'],
            m['enrichment'], m['pvalue'],
            m['mutations']))
    print()

file_list = ["FI_D1_D2_combined_WGS.maf.txt"]

if __name__ == '__main__':
    for filename in glob.glob(path):
    #for filename in file_list:
        vcf_files = [filename, ]
        with open("mut_and_sig.csv", "w") as csvfile:  # create csv file
            for file_name in vcf_files:
                motifs = get_motifs(file_name, assembly=37)
                print(file_name)
                print_motifs(motifs)
            filewriter = csv.writer(csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow([file_name])
            filewriter.writerow([motifs])
