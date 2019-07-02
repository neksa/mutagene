import requests
import json
import csv
import glob

MUTAGENE_URL = "https://www.ncbi.nlm.nih.gov/research/mutagene"

path = '/Users/cunninghamck/PycharmProjects/motifs/1st_analysis_tcga/*.maf.txt'
def get_profile(fname, assembly=37):
    """
    Calling MutaGene REST API to convert a VCF file into a mutational profile (96 context-dependent mutational probabilities)
    and profile_counts (counts of mutations for each of the 96 context-dependent mutations)
    It is important to specify genome assembly correctly. Curently 19, 37 and 38 will work
    """
    url = MUTAGENE_URL + '/pub/api/identify/profile'
    files = {'file': open(fname, 'rb')}
    r = requests.post(url, files=files, data={'assembly': assembly})
    # print("STATUS", r.status_code)
    if r.status_code == 200:
        return r.json()['result_counts']


def get_decomposition(profile_counts, signatures='COSMIC30'):
    """
    Decomposition of mutational profiles into a combination of signatures.
    It is highly recommended to use profile_counts instead of profile in order to use Maximum Likelihood method
    *signatures* should be one of COSMIC30  MUTAGENE5 MUTAGENE10
    *others_threshold* is used for not reporting signatures with exposure less or equal than the threshold and reporting the sum of their exposures as "Other signatures".
    Set *others_threshold* to 0 if not needed. The MutaGene website uses others_threshold = 0.05 by default.
    """
    url = MUTAGENE_URL + '/pub/api/identify/decomposition'
    r = requests.post(url, data={'profile_counts': json.dumps(profile_counts), 'signatures': signatures, 'others_threshold': 0.0})
    # print("STATUS", r.status_code)
    if r.status_code == 200:
        return r.json()['decomposition']


def print_profile_counts(profile_counts):
    """
    Printing context-dependent mutational profile
    """
    for mutation, value in profile.items():
        print("{}\t{:.0f}".format(mutation, value))
    print()


def print_decomposition(decomposition):
    """
    Printing the results of decomposition
    """
    for component in decomposition:
        print("{}\t{:.2f}\t{:.0f}".format(component['name'], component['score'], component['mutations']))
    print()


if __name__ == '__main__':
    for filename in glob.glob(path):
            vcf_files = [filename, ]
            with open("sig1.csv", "a") as csvfile:  # create csv file
                for file_name in vcf_files:
                    filewriter = csv.writer(csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL)
                    profile = get_profile(file_name, assembly=37)
                    print_profile_counts(profile)
                    if profile is not None:
                        for signature_type in ('COSMIC30', 'MUTAGENE5', 'MUTAGENE10'):
                            decomposition = get_decomposition(profile, signature_type)
                            print_decomposition(decomposition)
                    filewriter.writerow([file_name])
                    filewriter.writerow([decomposition])
                    filewriter.writerow([])
