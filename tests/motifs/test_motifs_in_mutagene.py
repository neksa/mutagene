
# from motifs_in_mutagene import identify_motifs, get_enrichment
# from mutations_io import read_mutations


# def test_motifs_in_mutagene():
#     with open("/panfs/pan1.be-md.ncbi.nlm.nih.gov/mutagene/data/skin/mutations/FI_D1_D2_combined_WGS.maf.txt") as f:
#         mut = f.read()
#         _, mutations_with_context, _ = read_mutations(mut, 'auto', 19)
#         matches = identify_motifs(mutations_with_context)
#         print(matches)


# def test_motifs_in_mutagene2():
#     mymotifs = [{
#         'name': 'Aristolochic Acid (AA)',
#         'logo': 'C[T>A]G',
#         'motif': 'T',
#         'position': 0,
#         'ref': 'T',
#         'alt': 'A',
#         'references': 'Hoang M. L. et al. Science Translational Medicine (2013)'
#     }]
#     with open("TCGA-2F-A9KP-01.maf.txt") as f:
#         mut = f.read()
#         _, mutations_with_context, _ = read_mutations(mut, 'auto', asm=37)
#         matches = identify_motifs(mutations_with_context, mymotifs)
#         print(matches)


# def test_motifs_in_mutagene3():
#     with open("/panfs/pan1.be-md.ncbi.nlm.nih.gov/mutagene/data/skin/mutations/FI_D1_D2_combined_WGS.maf.txt") as f:
#         _, mutations, _ = read_mutations(f.read(), 'MAF', 19)

#         observed2 = get_enrichment(mutations, "A", 0, "A", "T", 50)
#         print(observed2)


# if __name__ == '__main__':
#     #test_motifs_in_mutagene()
#    # print("-" * 60)
#     #print("-" * 60)
#     #test_motifs_in_mutagene2()
#     #print("-" * 60)
#     #print("-" * 60)
#     test_motifs_in_mutagene3()
#     #print("-" * 60)
#    # print("-" * 60)
